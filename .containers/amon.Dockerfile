ARG REF_NAME=main
FROM codecr.jlab.org/hallb/clas12/clas12-containers/hipo:main

RUN  apt install -y libgtkmm-4.0-dev
RUN echo ${REF_NAME} && git clone -b docker_build  https://code.jlab.org/touchte/amon.git \
    && cd  amon/src && make \
    && cp gui.exe /usr/local/bin/amon_gui





