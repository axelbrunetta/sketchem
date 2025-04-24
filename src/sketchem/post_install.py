# Run model installation if it hasn't run already when the package was installed

from .model_setup import install_decimer_model

def main():
    install_decimer_model()

if __name__ == "__main__":
    main()