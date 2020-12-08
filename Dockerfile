# EasyPQP Dockerfile
FROM python:3.9.1

# install numpy
RUN pip install numpy
RUN pip install pyprophet

# install EasyPQP and dependencies
ADD . /easypqp
WORKDIR /easypqp
RUN python setup.py install
WORKDIR /
RUN rm -rf /easypqp
