// Open image sequence in current directory from shell

cwd = call("java.lang.System.getProperty", "user.dir");
dir = getDirectory("current")
print("");
print("cwd:" + cwd );
print("dir:" + dir );

t = lengthOf(dir) == 0;
//print("t:" + t);

if (t)  {
   dir = cwd;
}
//print("dir:" + dir );

args = "open=[" + dir + "] file=tif sort use";
//print("args: " + args)

//print( "\n" );
//print( "image format: " + format );
//print( "path: " + impath );
//print( "image: " + getDirectory("image") );
//print( "temp: " + getDirectory("temp") );
//print( "home: " + getDirectory("home") );
//print( "startup: " + getDirectory("startup") );
//print( "imagej: " + getDirectory("imagej") );
//print("current: " + getDirectory("current"));

// Auto contrast stack
run("Appearance...", "no ij menu=0 gui=1 16-bit=Automatic");

// Open image stack
run("Image Sequence...", args);



