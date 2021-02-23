#!/bin/bash

set exuo pipefail

# some common dependencies ------------------------------------------------------------------------
sudo apt-get update
sudo apt-get install -y \
    build-essential curl git libfuse-dev fuse libcurl4-openssl-dev \
    curl wget software-properties-common libxml2-dev mime-support \
    automake libtool pkg-config libssl-dev ncurses-dev awscli \
    python-pip libbz2-dev liblzma-dev unzip imagemagick openjdk-11-jdk

# mount storage -----------------------------------------------------------------------------------
export TMPDIR=/mnt/local/temp
sudo mkfs -t ext4 /dev/nvme0n1 

sudo mkdir /mnt/local
sudo mount /dev/nvme0n1 /mnt/local
sudo chown ubuntu /mnt/local

mkdir /mnt/local/data
mkdir /mnt/local/temp
echo "export TMPDIR=/mnt/local/temp" >> ~/.profile

# conda setup -------------------------------------------------------------------------------------
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
eval "$($HOME/miniconda/bin/conda shell.bash hook)"
conda init
conda config --add channels bioconda

# tmux/neovim setup -------------------------------------------------------------------------------
echo "source-file ~/.tmux.d/.tmux.conf" > ~/.tmux.conf
git clone https://github.com/mchowdh200/.tmux.d.git ~/.tmux.d

sudo add-apt-repository ppa:neovim-ppa/stable -y
sudo apt-get update -y
sudo apt-get install -y neovim
git clone https://github.com/mchowdh200/.vim.git ~/.vim
mkdir ~/.config
mkdir ~/.config/nvim
printf "set runtimepath^=~/.vim runtimepath+=~/.vim/after\nlet &packpath=&runtimepath\nsource ~/.vim/vimrc" > ~/.config/nvim/init.vim
pip install jedi neovim
echo "alias vim=nvim" >> ~/.profile
echo "export EDITOR=nvim" >> ~/.profile

# setup path --------------------------------------------------------------------------------------
mkdir /mnt/local/bin
echo "PATH=$PATH:/mnt/local/bin" >> ~/.profile

# install gargs -----------------------------------------------------------------------------------
wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux -O /mnt/local/bin/gargs
chmod +x /mnt/local/bin/gargs

# install score client ----------------------------------------------------------------------------
wget -O score-client.tar.gz https://artifacts.oicr.on.ca/artifactory/dcc-release/bio/overture/score-client/[RELEASE]/score-client-[RELEASE]-dist.tar.gz
mkdir score-client &&
    tar -xvzf score-client.tar.gz -C score-client --strip-components 1
echo "PATH=$PATH:~/icgc-data/score-client/bin" >> ~/.profile

conda create -y -c bioconda -n smoove smoove

# install manta -----------------------------------------------------------------------------------
wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2

tar -xjvf manta-1.6.0.centos6_x86_64.tar.bz2 && 
    mv manta-1.6.0.centos6_x86_64 /mnt/local/manta/




