# Use an official Python runtime as a parent image.
FROM python:3.13.2-slim

# Install any system dependencies (e.g., gcc) that may be needed for building packages.
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
 && rm -rf /var/lib/apt/lists/*

# Set the working directory in the container.
WORKDIR /app

# Copy the requirements file into the container.
COPY requirements.txt .

# Install Python dependencies.
RUN pip install --upgrade pip && pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application code.
COPY . .

# Run the application.
CMD ["python", "calc_shadow.py"]