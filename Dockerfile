FROM python:3.8.2
EXPOSE 8501
WORKDIR ./
COPY requirement.txt ./requirement.txt
RUN pip3 install -r requirement.txt
COPY . .
CMD streamlit run chord_chart_2.py
