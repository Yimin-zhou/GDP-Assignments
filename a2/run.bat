@echo off
javac  -classpath jars\javaview.jar;jars\jvx.jar;jars\Jama-1.0.3.jar;. workshop\*.java
javac -classpath jars\javaview.jar;jars\jvx.jar;. menu\*.java
javac -classpath jars\javaview.jar;jars\jvx.jar;. test\*.java

@pause


set jv_jar=jars/javaview.jar;jars/jvx.jar;jars/vgpapp.jar;jars\Jama-1.0.3.jar;.
start java -cp %jv_jar% -Djava.library.path="dll" -Xmx1024m javaview model="models/bunny2000.obj" codebase=. archive.dev=show %*