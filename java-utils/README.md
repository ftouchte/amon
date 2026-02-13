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
  -Dfile=/home/touchte-codjo/Desktop/coatjava/coatjava/lib/clas/coat-libs-13.7.0-SNAPSHOT.jar \
  -DgroupId=org.jlab.coat \
  -DartifactId=coat-libs \
  -Dversion=13.7.0 \
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




