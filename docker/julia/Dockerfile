FROM julia:1.6
RUN apt-get -y update
RUN apt-get -y install git

RUN git clone https://github.com/cesaraustralia/GeneticRiskIndex /root/GeneticRiskIndex
RUN julia --project=/root/GeneticRiskIndex/julia -e "using Pkg; Pkg.instantiate()"

WORKDIR /root/GeneticRiskIndex/julia
USER root
CMD ["/bin/bash"]
