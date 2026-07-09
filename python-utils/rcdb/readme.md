# read_rcdb.py — Environnement Python isolé pour RCDB

Ce dépôt contient un script utilisant [RCDB](https://github.com/JeffersonLab/rcdb) (via l'installation
cvmfs de JLab, version `1.99.7`) pour interroger la base de données de runs.

## Problème rencontré

Le script plantait avec l'erreur suivante :

```
TypeError: Connection.execute() got an unexpected keyword argument 'run_max'
```

**Cause :** RCDB `1.99.7` a été écrit pour **SQLAlchemy 1.x**, où l'on peut passer les paramètres
de requête directement en `**kwargs` à `.execute()`. Dans **SQLAlchemy 2.x**, cette syntaxe n'est
plus supportée — les paramètres doivent être passés sous forme de dictionnaire. Comme
l'environnement Python système (cvmfs) installe SQLAlchemy 2.0.x, le code de RCDB casse.

**Solution :** utiliser un environnement virtuel Python (`venv`) isolé avec `sqlalchemy<2.0`,
sans toucher à l'installation système ni à d'autres projets.

## Prérequis

- Python 3 (testé avec Python 3.13)
- Accès en lecture au module RCDB installé sur cvmfs :
  `/u/scigroup/cvmfs/hallb/clas12/sw/noarch/rcdb/1.99.7/python/rcdb`

## Installation

```bash
# 1. Se placer dans le dossier du projet
cd /work/clas12/users/touchte/amon/python-utils

# 2. Créer un environnement virtuel isolé
python3 -m venv venv

# 3. Activer l'environnement
source venv/bin/activate

# 4. Installer les dépendances (versions compatibles avec RCDB 1.99.7)
pip install -r requirements.txt

# 5. Copier le module rcdb (cvmfs) directement dans le venv
#    -> évite tout problème de PYTHONPATH et toute interférence avec SQLAlchemy système
cp -r /u/scigroup/cvmfs/hallb/clas12/sw/noarch/rcdb/1.99.7/python/rcdb venv/lib/python3*/site-packages/
```

## Vérification

```bash
python -c "import sqlalchemy; print(sqlalchemy.__version__); print(sqlalchemy.__file__)"
python -c "import rcdb; print(rcdb.__file__)"
```

Les deux commandes doivent afficher un chemin contenant `venv/lib/...` — et **pas**
`/u/scigroup/cvmfs/...` — ce qui confirme que l'environnement est bien isolé du système.

`sqlalchemy.__version__` doit afficher une version `1.4.x`.

## Utilisation

À chaque nouvelle session de terminal :

```bash
cd /work/clas12/users/touchte/amon/python-utils
source venv/bin/activate
python read_rcdb.py
```

Pour quitter l'environnement virtuel :

```bash
deactivate
```

## Dépendances (`requirements.txt`)

```
sqlalchemy<2.0
pymysql
ply
```

## Notes

- Le dossier `venv/` est propre à chaque machine (chemins absolus, binaires compilés) et
  **ne doit pas être versionné sur Git** — il est exclu via `.gitignore`.
- Le module `rcdb` n'est pas installé via `pip` : il est copié depuis l'installation cvmfs
  de JLab directement dans le venv, à l'étape 5 ci-dessus. Si la version cvmfs de RCDB change
  de chemin, mettez à jour la commande `cp` en conséquence.
- Si d'autres erreurs `ModuleNotFoundError` apparaissent, installez les paquets manquants avec
  `pip install <nom_du_paquet>` (dans le venv activé), puis ajoutez-les à `requirements.txt`.