import ESRF_line as e
import click


@click.command()
@click.option("--n_part", default=1, help="Number of particle to simulate")
def main(n_part):
    e.run_simulation(n_part)


if __name__ == "__main__":
    main()
