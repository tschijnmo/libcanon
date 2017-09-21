# Dockerfile for a test environment for libcanon.

FROM gcc:latest

RUN apt-get update && apt-get install -y --no-install-recommends cmake

WORKDIR /libcanon/
ADD . .

