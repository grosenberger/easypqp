from datetime import datetime

import click


def timestamped_echo(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    click.echo(f"{timestamp} - {message}")
