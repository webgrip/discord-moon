# Use an official Python runtime as a parent image.
FROM python:3.14.5-slim@sha256:a7185a8e40af01bf891414a4df16ef10fc6000cee460a404a13da9029fe41604

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

COPY . .

# Run the application.
CMD ["python", "calc_shadow_spicey.py"]