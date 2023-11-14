FROM continuumio/miniconda3:4.12.0

### UPDATING CONDA ------------------------- ###

RUN conda update -y conda

### INSTALLING PIPELINE PACKAGES ----------- ###

# Installing CellRanger
# N.B. Download link expires fairly quickly, so has to be changed
#RUN cd /opt && \
#    curl -o cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.1.tar.gz?Expires=1664073344&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMS50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NjQwNzMzNDR9fX1dfQ__&Signature=lPqgTcmtYTJAxlYO-76O4vDzQ7hPp4as8SZOHfaJEmWtKDYMGDTaC4Gd37~N~0pvBqWkxyfdvaj4HCIk0svpk3qLaIqBaTkSxNbs~jag6mo6CVKqq-D8~0ulbkONly2EdbAnv8zEM8LTjA6kaWuUpR5SC~bDSnGhr2N2RrT7U-N-xjyzxRViX2ONZbe5KkoX4HHbm5n0a3yZwzmjmSwgxpB1S9jtb~3U2T~GpfBTtt19M82Oif5lwNvbjL-TIcZ44jVdBJE2q~ix22unpUNuSlfh~nXUzv6TcfaQWrLdt5~ssUuEo1YmS6QR5rLY7jOFDvo~~ySXpwEpb0rlTWoCzg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \
#    tar -xzvf cellranger-7.1.0.tar.gz && \
#    rm -f cellranger-7.1.0.tar.gz

# Alternatively, if a copy of CellRanger was already downloaded, use the following
COPY cellranger-7.1.0 /opt/cellranger-7.1.0

# Add Cellranger to path
ENV PATH $PATH:/opt/cellranger-7.1.0

### SETTING WORKING ENVIRONMENT ------------ ###

# Set workdir to /home/
WORKDIR /home/

# Launch bash automatically
CMD ["/bin/bash"]
