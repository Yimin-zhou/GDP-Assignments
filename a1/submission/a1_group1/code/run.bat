@echo off
set jv_jar=jars/javaview.jar;jars/jvx.jar;jars/vgpapp.jar;jars\Jama-1.0.3.jar;.
start java -cp %jv_jar% -Djava.library.path="dll" -Xmx1024m javaview model="models/caesarRange_4.obj.zip" codebase=. archive.dev=show %*