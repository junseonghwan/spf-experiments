#gsutil -m cp build/distributions/spf-experiments.zip gs://seonghwanjun_spf/
scp -i "spf.pem" build/distributions/spf-experiments.zip ubuntu@$1:
