import click
from gevent.pywsgi import WSGIServer
from gevent import monkey
from mmsplice.api import app


@click.group()
def cli():
    pass


@cli.command()
@click.option('--port', default=5000, help='Port of python server')
def run_api(port):
    monkey.patch_all()
    http_server = WSGIServer(('', port), app.wsgi_app)
    http_server.serve_forever()


if __name__ == '__main__':
    cli()
