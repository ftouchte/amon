# Dependencies

One needs to install coatjava and `ahdc/alignment-v?` branch

```shell
get clone https://github.com/ftouchte/coatjava.git
cd coatjava
git checkout ahdc/alignment-v?
```

Each time the modifications in **jeffersonlab/coatjava:development** will be merge in the `ahdc/alignment-v?` branch, one needs to config the path in the `pom.xml` file

1. Look at the `pom.xml` file
1. Identify all references to local jar files with `coatjava` in it
1. Look at the corresponding directory paths and use the file version that is actually installed (it generally results in just changing the version number, e.g moving from 13.0.0 to 14.1.1)

```shell
# In vim 
vim pomw.xml
:%s/13.0.0/14.1.1/g

# After that, in the terminal, in the directory where is located this pom.xml file
mvn install
```


# Create a project with maven

```shell
mvn archetype:generate \
  -DgroupId=io.github.ftouchte \
  -DartifactId=java-utils \
  -DarchetypeArtifactId=maven-archetype-quickstart \
  -DinteractiveMode=false
```

# Add an external dependency in the pom.xml

## If a maven repository exists

```xml
<!-- Source: https://mvnrepository.com/artifact/org.jlab.coat/coat-libs -->
<dependency>
    <groupId>org.jlab.coat</groupId>
    <artifactId>coat-libs</artifactId>
    <version>13.0.0</version>
    <scope>compile</scope>
</dependency>
```

## Local installation from a jar file

Example

```shell
mvn install:install-file \
  -Dfile=/w/hallb-scshelf2102/clas12/users/touchte/coatjava/coatjava/lib/clas/coat-libs-13.7.1-SNAPSHOT.jar \
  -DgroupId=org.jlab.coat \
  -DartifactId=coat-libs \
  -Dversion=13.7.1 \
  -Dpackaging=jar
```

add the dependency in the file `pom.xml`

```xml
<dependency>
    <groupId>org.jlab.coat</groupId>
    <artifactId>coat-libs</artifactId>
    <version>13.7.0</version>
</dependency>
```

## Find dependencies

```
# Find the jar files of a specific dependency. Example:
mvn dependency:build-classpath | grep "xchart"
```





