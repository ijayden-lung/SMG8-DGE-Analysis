# SMG8-DGE-Analysis


## Working environment

### Use `jupyter/datascience-notebook`

Plan to use 
`docker pull jupyter/datascience-notebook`([link](https://hub.docker.com/r/jupyter/datascience-notebook)) as a base working environment. It has combined two docker images: `jupyter/scipy-notebook` and `jupyter/r-notebook`. A `jupyter/all-spark-notebook` in "[Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html)" is also worth to check out if current base is not enough. (I actually checked `jupyter/all-spark-notebook`, `R sessionInfo()` is the same as that in `jupyter/datascience-notebook`. The python in them show slight differences like `spark` related. One major difference is in `all-spark-notebook` "numpy==1.16.1", while in `datascience-notebook` it is "numpy==1.13.3". I decided to use `datascience-notebook` first, as the name indicates "data science", which supposed to have been runed to fit data science.)

### Launch it for the first time

Create a data folder on host and modify its permission ([reason about permission](https://github.com/jupyter/docker-stacks/issues/114#issuecomment-242682549)). This is only needed for the first time. 

```
mkdir jupyterLab
chmod 777 jupyterLab
```

Pull the docker image and launch it. This is only needed for the first time. 

```
docker pull jupyter/datascience-notebook
docker run -d -P -v /home/data/jupyterLab:/home/jovyan/work --name jupyter_datascience jupyter/datascience-notebook
```

### Stop and restart it
```
docker stop jupyter_datascience # stop the instance
docker start jupyter_datascience
docker port jupyter_datascience # find the correct port number on host
docker logs --tail 3 jupyter_datascience # find the token to login
```

Put all notebooks under `work` to get them saved on host, so that once you restart the instance, they are still there.

### Notebooks

All the notebooks for analysis are under [notebooks/](https://github.com/RodenLuo/RBP-PAS/tree/master/notebooks)

# Reports

All the reports are under [reports](https://github.com/RodenLuo/RBP-PAS/tree/master/reports)
