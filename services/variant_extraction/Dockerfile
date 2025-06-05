FROM python:3.11-slim
WORKDIR /app
COPY app/ app/
COPY requirements.txt .
RUN pip install -r requirements.txt
CMD ["python", "app/parser.py"]
