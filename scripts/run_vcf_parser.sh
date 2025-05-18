#!/bin/bash
docker-compose build vcf-parser
docker-compose run --rm vcf-parser
