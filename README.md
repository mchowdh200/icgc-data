# Run setup script
```
bash icgc_s3_score_client_setup.sh
```

# Add your access token
edit `icgc-data/score-client/conf/application.properties`
```bash
accessToken= # add your token here
```

# Untar the manifest
```bash
tar -xvzf manifest.1599710870723.tar.gz
```

# Download bams using score-client and manifest
```
score-client/bin/score-client download --manifest manifest.aws-virginia.1599710870723.tsv --output-dir /mnt/local/data
```
