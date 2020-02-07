# EasyPQP Dockerfile
FROM python:3.7.3

# install numpy
RUN pip install numpy pyprophet

# install EasyPQP and dependencies
ADD . /easypqp
WORKDIR /easypqp
RUN python setup.py install
WORKDIR /
RUN rm -rf /easypqp
