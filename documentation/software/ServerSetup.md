Servers were run on Ubuntu 24.04

After cloning the webserver, run the following commands to obtain a docker container


Installing docker:
```
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt update

sudo apt install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

After installing docker, before running it, we need to create an `.env.local` file, that would be composed of these values 

```
BREVO_API_KEY=... # Key for the BREVO API that is used to send emails
FROM_EMAIL=... # Email you have used for BREVO
FROM_NAME=... # Name that would appear in the email
```