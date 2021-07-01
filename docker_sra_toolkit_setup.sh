#!/bin/env bash

# install docker
sudo apt-get install apt-transport-https ca-certificates curl gnupg lsb-release
curl -fsSL https://download.docker.com/linux/ubuntu/gpg |
    sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo "deb [arch=amd64 signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] \
 https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" |
    sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io

# move the docker location before the stuff below
sed 's/ExecStart=\/usr\/bin\/dockerd -H fd:\/\/ --containerd=\/run\/containerd\/containerd.sock/ExecStart=\/usr\/bin\/dockerd -H fd:\/\/ --containerd=\/run\/containerd\/containerd.sock -g \/mnt\/local\/docker/g' \
    /lib/systemd/system/docker.service > temp.txt
sudo mv temp.txt /lib/systemd/system/docker.service
sudo systemctl stop docker
sudo systemctl daemon-reload
mkdir /mnt/local/docker
sudo systemctl start docker
sudo usermod -aG docker ubuntu
newgrp docker

docker pull ncbi/sra-tools
