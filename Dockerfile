FROM allenlawrence94/jupyter-base:latest

#RUN conda config --set auto_update_conda false
RUN conda install -c conda-forge cvxopt --yes

RUN conda install pandas numpy --yes
