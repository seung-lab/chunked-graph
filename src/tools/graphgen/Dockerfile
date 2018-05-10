FROM python:3-slim-stretch
LABEL MAINTAINER="Nico Kemnitz"

ENV JULIA_PATH /usr/local/julia
ENV JULIA_VERSION 0.6.2

RUN set -ex \
	&& BUILD_PACKAGES=" \
		curl \
		gcc \
		g++ \
		git \
		make \
	" \
	&& apt-get update && apt-get install -y $BUILD_PACKAGES --no-install-recommends && rm -rf /var/lib/apt/lists/* \
	\
	&& curl -fL -o julia.tar.gz "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VERSION%[.-]*}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" \
	&& mkdir "$JULIA_PATH" \
	&& tar -xzf julia.tar.gz -C "$JULIA_PATH" --strip-components 1 \
	&& rm julia.tar.gz \
	&& ln -s "$JULIA_PATH"/bin/julia /usr/bin/julia \
	\
	&& pip install --no-cache-dir cloud-volume \
	&& mkdir /root/.cloudvolume \
	&& echo "" >> /root/.cloudvolume/project_name \
	&& ln -s /secrets /root/.cloudvolume/secrets \
	\
	&& julia -e "ENV[\"PYTHON\"] = \"python\"; Pkg.add(\"PyCall\");" \
	&& julia -e "Pkg.clone(\"https://github.com/seung-lab/CloudVolume.jl\"); using CloudVolume" \
	&& julia -e "Pkg.clone(\"https://github.com/nkemnitz/MySQL.jl\"); Pkg.checkout(\"MySQL\", \"binary-length-check\"); Pkg.build(\"MySQL\"); using MySQL" \
	\
	&& git clone "https://github.com/seung-lab/chunked-graph.git" /usr/local/chunked-graph \
	\
	&& julia -e "Pkg.add(\"DataStructures\"); using DataStructures" \
	&& julia -e "Pkg.add(\"IterTools\"); using IterTools" \
	&& julia -e "Pkg.add(\"LightGraphs\"); using LightGraphs" \
	&& julia -e "Pkg.add(\"LightGraphsFlows\"); using LightGraphsFlows" \
	&& julia -e "Pkg.add(\"Memento\"); using Memento" \
	&& julia -e "Pkg.add(\"CodecZstd\"); using CodecZstd" \
	&& julia -e "Pkg.add(\"DataFrames\"); using DataFrames" \
	&& julia -e "Pkg.add(\"Memoize\"); using Memoize" \
	&& julia -e "Pkg.add(\"NamedTuples\"); using NamedTuples" \
	&& echo "Skipping julia -e \"Pkg.add(\"MySQL\"); using MySQL\"" \
	&& julia -e "Pkg.add(\"OffsetArrays\"); using OffsetArrays" \
	\
	&& apt-get purge -y --auto-remove $BUILD_PACKAGES \
	&& find /usr/local -depth \
		\( \
			\( -type d -a \( -name __pycache__ -o -name test -o -name tests \) \) \
			-o \
			\( -type f -a \( -name '*.pyc' -o -name '*.pyo' \) \) \
		\) -exec rm -rf '{}' + \
	&& rm -rf /usr/src/python \
	\
	&& find /root/.julia/v0.6 -depth \
		\( \
			\( -type d -a \( -name test -o -name tests -o -name .git \) \) \
			-o \
			\( -type f -a \( -name '*.gitignore' -o -name '*.jl.cov' -o -name '*.jl.mem' \) \) \
		\) -exec rm -rf '{}' + \
	&& rm -rf /root/.julia/v0.6/METADATA \
	&& rm -rf /root/.julia/.cache \
	&& rm -rf /usr/local/julia/share/doc \
	\
	&& rm -rf /var/lib/apt/lists/*
