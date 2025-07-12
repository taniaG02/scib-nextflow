#!/bin/bash

# Activar Conda
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
else
    echo "No se encontró la instalación de Miniconda en $HOME/miniconda3"
    exit 1
fi

# Carpeta donde están los YAML
ENV_DIR="conda"

# Crear entornos desde cada archivo .yml
for yml in "$ENV_DIR"/*.yml; do
    env_name=$(basename "$yml" .yml)
    echo "Instalando entorno: $env_name desde $yml"
    conda env create -f "$yml" -n "$env_name"
done

echo "Todos los entornos han sido instalados localmente."
