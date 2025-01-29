# Installation

```{note}
This guide explains how to install **dem_lic**, including dependencies and post-installation checks.
```

## Prerequisites

Before installing **dem_lic**, ensure you have:

- **Python 3.8+** installed. Check your version with:

  ```bash
  python --version
  ```

- **pip** and **virtual environment tools** installed:

  ```bash
  python -m ensurepip --default-pip
  python -m pip install --upgrade pip virtualenv
  ```


## Installation Steps

### Clone the Repository

You can install **dem_lic** by cloning the GitHub repository.

#### HTTP Method

```bash
git clone https://github.com/ESaint-Denis/dem_lic.git
cd dem_lic
pip install -e .
```

#### SSH Method

```bash
git clone git@github.com:ESaint-Denis/dem_lic.git
cd dem_lic
pip install -e .
```


## Post-Installation Verification

To confirm the package is correctly installed, run:

```bash
dem_lic --help
```

