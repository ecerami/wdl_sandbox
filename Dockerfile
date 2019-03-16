FROM python:2.7

RUN pip install pytest
RUN apt-get update && apt-get install -y curl && apt-get clean && rm -rf /var/lib/apt/lists

