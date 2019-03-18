FROM ctzhu/conda3-mipcl

RUN conda install jupyter --yes

RUN mkdir /root/.jupyter && \
    mkdir /root/work

ADD jupyter_notebook_config.py /root/.jupyter/

WORKDIR /root/work

CMD jupyter notebook

#RUN conda config --set auto_update_conda false
RUN conda install -c conda-forge cvxopt --yes

RUN conda install pandas numpy matplotlib --yes
