# import tensorflow as tf
import numpy as np
import tensorflow.keras.backend as K
from tensorflow.keras.layers import Layer
from tensorflow.keras.layers import Conv1D
from tensorflow.keras.regularizers import Regularizer
from tensorflow.keras import initializers
import scipy.interpolate as si


DNA = ["A", "C", "G", "T"]


def normalize_data_format(value):
    if value is None:
        value = K.image_data_format()
    data_format = value.lower()
    if data_format not in {'channels_first', 'channels_last'}:
        raise ValueError('The `data_format` argument must be one of '
                         '"channels_first", "channels_last". Received: ' +
                         str(value))
    return data_format


class GlobalAveragePooling1D_Mask0(Layer):
    """
    Global average pooling operation for temporal data.
    Masking out 0-padded input.
    """

    def __init__(self, data_format='channels_last', **kwargs):
        super(GlobalAveragePooling1D_Mask0, self).__init__(**kwargs)
        self.data_format = normalize_data_format(data_format)

    def compute_output_shape(self, input_shape):
        input_shape = input_shape[0]
        if self.data_format == 'channels_first':
            return (input_shape[0], input_shape[1])
        else:
            return (input_shape[0], input_shape[2])

    def call(self, inputs):
        inputs, model_inputs = inputs
        steps_axis = 1 if self.data_format == 'channels_last' else 2
        mask = K.max(model_inputs, axis=2, keepdims=True)
        inputs *= mask
        return K.sum(inputs, axis=steps_axis) / K.maximum(
            K.sum(mask, axis=steps_axis), K.epsilon())


class ConvSequence(Conv1D):
    VOCAB = DNA

    def __init__(self,
                 filters,
                 kernel_size,
                 strides=1,
                 padding='valid',
                 dilation_rate=1,
                 activation=None,
                 use_bias=True,
                 kernel_initializer='glorot_uniform',
                 bias_initializer='zeros',
                 kernel_regularizer=None,
                 bias_regularizer=None,
                 activity_regularizer=None,
                 kernel_constraint=None,
                 bias_constraint=None,
                 seq_length=None,
                 **kwargs):

        # override input shape
        if seq_length:
            kwargs["input_shape"] = (seq_length, len(self.VOCAB))
            kwargs.pop("batch_input_shape", None)

        super(ConvSequence, self).__init__(
            filters=filters,
            kernel_size=kernel_size,
            strides=strides,
            padding=padding,
            dilation_rate=dilation_rate,
            activation=activation,
            use_bias=use_bias,
            kernel_initializer=kernel_initializer,
            bias_initializer=bias_initializer,
            kernel_regularizer=kernel_regularizer,
            bias_regularizer=bias_regularizer,
            activity_regularizer=activity_regularizer,
            kernel_constraint=kernel_constraint,
            bias_constraint=bias_constraint,
            **kwargs)

        self.seq_length = seq_length

    def build(self, input_shape):
        if int(input_shape[-1]) != len(self.VOCAB):
            raise ValueError("{cls} requires input_shape[-1] == {n}. Given: {s}".
                             format(cls=self.__class__.__name__, n=len(self.VOCAB), s=input_shape[-1]))
        return super(ConvSequence, self).build(input_shape)

    def get_config(self):
        config = super(ConvSequence, self).get_config()
        config["seq_length"] = self.seq_length
        return config


class ConvDNA(ConvSequence):
    VOCAB = DNA
    VOCAB_name = "DNA"


def get_S(n_bases=10, spline_order=3, add_intercept=True):
    # mvcv R-code
    # S<-diag(object$bs.dim);
    # if (m[2]) for (i in 1:m[2]) S <- diff(S)
    # object$S <- list(t(S)%*%S)  # get penalty
    # object$S[[1]] <- (object$S[[1]]+t(object$S[[1]]))/2 # exact symmetry

    S = np.identity(n_bases)
    m2 = spline_order - 1  # m[2] is the same as m[1] by default

    # m2 order differences
    for i in range(m2):
        S = np.diff(S, axis=0)  # same as diff() in R

    S = np.dot(S.T, S)
    S = (S + S.T) / 2  # exact symmetry

    if add_intercept is True:
        # S <- cbind(0, rbind(0, S)) # in R
        zeros = np.zeros_like(S[:1, :])
        S = np.vstack([zeros, S])

        zeros = np.zeros_like(S[:, :1])
        S = np.hstack([zeros, S])

    return S.astype(np.float32)


def get_knots(start, end, n_bases=10, spline_order=3):
    """
    Arguments:
        x; np.array of dim 1
    """
    x_range = end - start
    start = start - x_range * 0.001
    end = end + x_range * 0.001

    # mgcv annotation
    m = spline_order - 1
    nk = n_bases - m            # number of interior knots

    dknots = (end - start) / (nk - 1)
    knots = np.linspace(start=start - dknots * (m + 1),
                        stop=end + dknots * (m + 1),
                        num=nk + 2 * m + 2)
    return knots.astype(np.float32)


def get_X_spline(x, knots, n_bases=10, spline_order=3, add_intercept=True):
    """
    Returns:
        np.array of shape [len(x), n_bases + (add_intercept)]
    # BSpline formula
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html#scipy.interpolate.BSpline
    Fortran code:
    https://github.com/scipy/scipy/blob/v0.19.0/scipy/interpolate/fitpack/splev.f
    """
    if len(x.shape) is not 1:
        raise ValueError("x has to be 1 dimentional")

    tck = [knots, np.zeros(n_bases), spline_order]

    X = np.zeros([len(x), n_bases])

    for i in range(n_bases):
        vec = np.zeros(n_bases)
        vec[i] = 1.0
        tck[1] = vec

        X[:, i] = si.splev(x, tck, der=0)

    if add_intercept is True:
        ones = np.ones_like(X[:, :1])
        X = np.hstack([ones, X])

    return X.astype(np.float32)


class BSpline():
    """Class for computing the B-spline funcions b_i(x) and
    constructing the penality matrix S.
    # Arguments
        start: float or int; start of the region
        end: float or int; end of the region
        n_bases: int; number of spline bases
        spline_order: int; spline order
    # Methods
        - **getS(add_intercept=False)** - Get the penalty matrix S
              - Arguments
                     - **add_intercept**: bool. If true, intercept column is added to the returned matrix.
              - Returns
                     - `np.array`, of shape `(n_bases + add_intercept, n_bases + add_intercept)`
        - **predict(x, add_intercept=False)** - For some x, predict the bn(x) for each base
              - Arguments
                     - **x**: np.array; Vector of dimension 1
                     - **add_intercept**: bool; If True, intercept column is added to the to the final array
              - Returns
                     - `np.array`, of shape `(len(x), n_bases + (add_intercept))`
    """

    def __init__(self, start=0, end=1, n_bases=10, spline_order=3):

        self.start = start
        self.end = end
        self.n_bases = n_bases
        self.spline_order = spline_order

        self.knots = get_knots(self.start, self.end,
                               self.n_bases, self.spline_order)

        self.S = get_S(self.n_bases, self.spline_order, add_intercept=False)

    def __repr__(self):
        return "BSpline(start={0}, end={1}, n_bases={2}, spline_order={3})".\
            format(self.start, self.end, self.n_bases, self.spline_order)

    def getS(self, add_intercept=False):
        """Get the penalty matrix S
        Returns
            np.array, of shape (n_bases + add_intercept, n_bases + add_intercept)
        """
        S = self.S
        if add_intercept is True:
            # S <- cbind(0, rbind(0, S)) # in R
            zeros = np.zeros_like(S[:1, :])
            S = np.vstack([zeros, S])

            zeros = np.zeros_like(S[:, :1])
            S = np.hstack([zeros, S])
        return S

    def predict(self, x, add_intercept=False):
        """For some x, predict the bn(x) for each base
        Arguments:
            x: np.array; Vector of dimension 1
            add_intercept: bool; should we add the intercept to the final array
        Returns:
            np.array, of shape (len(x), n_bases + (add_intercept))
        """
        # sanity check
        if x.min() < self.start:
            raise Warning("x.min() < self.start")
        if x.max() > self.end:
            raise Warning("x.max() > self.end")

        return get_X_spline(x=x,
                            knots=self.knots,
                            n_bases=self.n_bases,
                            spline_order=self.spline_order,
                            add_intercept=add_intercept)

    def get_config(self):
        return {"start": self.start,
                "end": self.end,
                "n_bases": self.n_bases,
                "spline_order": self.spline_order
                }

    @classmethod
    def from_config(cls, config):
        return cls(**config)


class GAMRegularizer(Regularizer):

    def __init__(self, n_bases=10, spline_order=3, l2_smooth=0., l2=0.):
        """Regularizer for GAM's
        # Arguments
            n_bases: number of b-spline bases
            order: spline order (2 for quadratic, 3 for qubic splines)
            l2_smooth: float; Smoothness penalty (penalize w' * S * w)
            l2: float; L2 regularization factor - overall weights regularizer
        """
        # convert S to numpy-array if it's a list

        self.n_bases = n_bases
        self.spline_order = spline_order
        self.l2_smooth = K.cast_to_floatx(l2_smooth)
        self.l2 = K.cast_to_floatx(l2)

        # convert to K.constant
        self.S = K.constant(
            K.cast_to_floatx(
                get_S(n_bases, spline_order, add_intercept=False)
            ))

    def __call__(self, x):
        # x.shape = (n_bases, n_spline_tracks)
        # from conv: (kernel_width=1, n_bases, n_spline_tracks)
        from_conv = len(K.int_shape(x)) == 3
        if from_conv:
            x = K.squeeze(x, 0)

        n_spline_tracks = K.cast_to_floatx(K.int_shape(x)[1])

        regularization = 0.

        if self.l2:
            regularization += K.sum(self.l2 * K.square(x)) / n_spline_tracks

        if self.l2_smooth:
            # https://keras.io/backend/#batch_dot
            # equivalent to mean( diag(x' * S * x) )
            regularization += self.l2_smooth * \
                K.mean(K.batch_dot(x, K.dot(self.S, x), axes=1))

        return regularization

    def get_config(self):
        # convert S to list()
        return {'n_bases': self.n_bases,
                'spline_order': self.spline_order,
                'l2_smooth': float(self.l2_smooth),
                'l2': float(self.l2),
                }


class SplineWeight1D(Layer):
    """Up- or down-weight positions in the activation array of 1D convolutions:
    `x^{out}_{ijk} = x^{in}_{ijk}* (1 + f_S^k(j)) \;,`
    where f_S is the spline transformation.
    # Arguments
        n_bases: int; Number of spline bases used for the positional effect.
        l2_smooth: (float) L2 regularization strength for the second
    order differences in positional bias' smooth splines. (GAM smoothing regularization)
        l2: (float) L2 regularization strength for the spline base coefficients.
        use_bias: boolean; should we add a bias to the transition
        bias_initializer: bias initializer - from `keras.initializers`
    """

    def __name__(self):
        return "SplineWeight1D"

    def __init__(self,
                 # spline type
                 n_bases=10,
                 spline_degree=3,
                 share_splines=False,
                 # regularization
                 l2_smooth=0,
                 l2=0,
                 use_bias=False,
                 bias_initializer='zeros',
                 **kwargs):
        self.n_bases = n_bases
        self.spline_degree = spline_degree
        self.share_splines = share_splines
        self.l2 = l2
        self.l2_smooth = l2_smooth
        self.use_bias = use_bias
        self.bias_initializer = initializers.get(bias_initializer)

        super(SplineWeight1D, self).__init__(**kwargs)

    def build(self, input_shape):
        # input_shape = (None, steps, filters)

        start = 0
        end = int(input_shape[1])
        filters = int(input_shape[2])

        if self.share_splines:
            n_spline_tracks = 1
        else:
            n_spline_tracks = filters

        # setup the bspline object
        self.bs = BSpline(start, end - 1,
                          n_bases=self.n_bases,
                          spline_order=self.spline_degree
                          )

        # create X_spline,
        self.positions = np.arange(end)
        # shape = (end, self.n_bases)
        self.X_spline = self.bs.predict(self.positions, add_intercept=False)

        # convert to the right precision and K.constant
        self.X_spline_K = K.constant(K.cast_to_floatx(self.X_spline))

        # add weights - all set to 0
        self.kernel = self.add_weight(shape=(self.n_bases, n_spline_tracks),
                                      initializer='zeros',
                                      name='kernel',
                                      regularizer=GAMRegularizer(self.n_bases, self.spline_degree,
                                                                 self.l2_smooth, self.l2),
                                      trainable=True)

        if self.use_bias:
            self.bias = self.add_weight((n_spline_tracks, ),
                                        initializer=self.bias_initializer,
                                        name='bias',
                                        regularizer=None)

        # Be sure to call this somewhere!
        super(SplineWeight1D, self).build(input_shape)

    def call(self, x):

        spline_track = K.dot(self.X_spline_K, self.kernel)

        if self.use_bias:
            spline_track = K.bias_add(spline_track, self.bias)

        # if self.spline_exp:
        #     spline_track = K.exp(spline_track)
        # else:
        spline_track = spline_track + 1

        # multiply together the two coefficients
        output = spline_track * x

        return output

    def compute_output_shape(self, input_shape):
        return input_shape

    def get_config(self):
        config = {
            'n_bases': self.n_bases,
            'spline_degree': self.spline_degree,
            'share_splines': self.share_splines,
            # 'spline_exp': self.spline_exp,
            'l2_smooth': self.l2_smooth,
            'l2': self.l2,
            'use_bias': self.use_bias,
            'bias_initializer': initializers.serialize(self.bias_initializer),
        }
        base_config = super(SplineWeight1D, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))

    def positional_effect(self):
        w = self.get_weights()[0]
        pos_effect = np.dot(self.X_spline, w)
        return {"positional_effect": pos_effect, "positions": self.positions}
