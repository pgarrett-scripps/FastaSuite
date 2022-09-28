FROM python:3.10.7-alpine3.16

RUN apk add git

WORKDIR /usr/src/app

COPY requirements.txt ./

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD streamlit run ./home.py --server.maxUploadSize 2000