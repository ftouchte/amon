# Utiliser Ubuntu comme image de base
FROM ubuntu:20.04

# Installer les outils nécessaires
RUN apt-get update && \
    apt-get install -y g++ cmake make && \
    rm -rf /var/lib/apt/lists/*

# Définir le répertoire de travail
WORKDIR /app

# Copier les fichiers du projet
COPY . .

# Compiler le code source
RUN g++ -o mon_application src/main.cpp

# Commande par défaut lors de l'exécution du conteneur
CMD ["./mon_application"]
