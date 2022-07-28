#!/bin/bash

set exuo pipefail


### some common dependencies
################################################################################
sudo apt-get update
sudo apt-get install -y \
    build-essential curl git libcurl4-openssl-dev \
    curl wget software-properties-common libxml2-dev mime-support \
    automake libtool pkg-config libssl-dev ncurses-dev awscli \
    python-pip libbz2-dev liblzma-dev unzip openjdk-11-jdk

### mount storage
################################################################################
export TMPDIR=/mnt/local/temp
sudo mkdir /mnt/local

# setup raid 0 if more than one drive specified
# the nvme drive naming convention is not consistent enough
# so I have just resorted to filtering out nvme disks with
# the expected size (ex smallest c5d has a 50GB so > 40 is my threshold)
num_drives=$(lsblk -b -o NAME,SIZE | grep 'nvme'| awk '$2 > 4e10' | wc -l)
drive_list=$(lsblk -b -o NAME,SIZE | grep 'nvme' |
                awk '$2 > 4e10' |
                awk 'BEGIN{ORS=" "}{print "/dev/"$1 }')
if [[ $num_drives > 1 ]]; then
    drive_list=$(lsblk -b -o NAME,SIZE | grep 'nvme' |
                 awk '$2 > 4e10' |
                 awk 'BEGIN{ORS=" "}{print "/dev/"$1 }')
    sudo mdadm --create --verbose \
         /dev/md0 \
         --level=0 \
         --raid-devices=$num_drives $drive_list
    sudo mkfs -t xfs /dev/md0
    sudo mount /dev/md0 /mnt/local
else
    sudo mkfs -t xfs $drive_list
    sudo mount $drive_list /mnt/local
fi

user=$(whoami)
sudo chown $user /mnt/local
mkdir /mnt/local/data
mkdir /mnt/local/temp
echo "export TMPDIR=/mnt/local/temp" >> ~/.profile

### conda/mamba setup
################################################################################
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /mnt/local/miniconda.sh
bash /mnt/local/miniconda.sh -b -p /mnt/local/miniconda
eval "$(/mnt/local/miniconda/bin/conda shell.bash hook)"
conda init
conda config --add channels bioconda
conda install -y -c conda-forge mamba


### tmux/neovim setup
################################################################################
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

### setup path/other environment variables
################################################################################
mkdir /mnt/local/bin
echo 'PATH=$PATH:/mnt/local/bin' >> ~/.profile
echo "PS1='\[\e[01;32m\]\u@\h\[\e[0m\] \[\e[34m\]\w\[\e[0m\] \n$ '" >> ~/.profile

### install gargs
################################################################################
wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux -O /mnt/local/bin/gargs
chmod +x /mnt/local/bin/gargs

### install score client
################################################################################
wget -O score-client.tar.gz https://artifacts.oicr.on.ca/artifactory/dcc-release/bio/overture/score-client/[RELEASE]/score-client-[RELEASE]-dist.tar.gz
mkdir score-client &&
    tar -xvzf score-client.tar.gz -C score-client --strip-components 1
echo 'PATH=$PATH:~/icgc-data/score-client/bin' >> ~/.profile
