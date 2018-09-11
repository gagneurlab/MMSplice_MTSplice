import click
from mmsplice.api import app


@click.group()
def cli():
    pass


@cli.command()
def run_api():
    app.run(debug=True)


if __name__ == '__main__':
    cli()
