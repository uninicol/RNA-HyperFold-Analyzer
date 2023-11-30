FROM python:3.10.13

WORKDIR /rna-folding-hypergraph

COPY requirements.txt ./requirements.txt

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["python3", "src/main.py"]