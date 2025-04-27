![Project Logo](assets/banner.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
sketchem
</h1>

<br>


A molecular pictionary



## 👩‍💻 Installation

Create a new environment, you may also give the environment a different name.

```
conda create -n sketchem python=3.10
```

```
conda activate sketchem
(conda_env) $ pip install .
```

The DECIMER Canonical model will be automatically downloaded during installation. If you need to manually trigger the model download, you can run:

```bash
sketchem-setup
```


If you need jupyter lab, install it 

```
(sketchem) $ pip install jupyterlab
```

## How to deploy on Streamlit Cloud


### Setup Streamlit Cloud

1. Fork or push this repository to your GitHub account
2. Connect your repository to Streamlit Cloud
3. Add the following secret in the Streamlit Cloud dashboard:

```
GEMINI_API_KEY = "your-gemini-api-key"
```

### Local Development

keep in mind that mulptiplayer game mode is not available when running the app locally, it only works on the deployed version


1. Create a `.env` file in the root of your project with:
```
GEMINI_API_KEY=your-gemini-api-key
```

2. Run the app:
```
streamlit run src/sketchem/main.py
```

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```



