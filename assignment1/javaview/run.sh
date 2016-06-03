#!/bin/bash
jv_jar=jars/javaview.jar:jars/jvx.jar:jars/vgpapp.jar:jars\Jama-1.0.3.jar:commons-math3-3.6.1/commons-math3-3.6.1.jar:.
java -cp $jv_jar -Djava.library.path="dll" -Xmx1024m javaview model='models/rabbit_simple.obj' codebase=. #archive.dev=show %*
