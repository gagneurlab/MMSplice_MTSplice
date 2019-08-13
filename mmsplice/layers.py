from keras.layers.pooling import Layer
import keras.backend as K


class GlobalAveragePooling1D_Mask0(Layer):
    """
    Global average pooling operation for temporal data.
    Masking out 0-padded input.
    """

    def __init__(self, data_format='channels_last', **kwargs):
        super(GlobalAveragePooling1D_Mask0, self).__init__(**kwargs)
        self.data_format = K.normalize_data_format(data_format)

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
