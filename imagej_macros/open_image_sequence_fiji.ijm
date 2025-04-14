// Open image sequence in current directory from shell

cwd = call("java.lang.System.getProperty", "user.dir");
cur = getDirectory("current")
//print("");
//print("cwd:" + cwd );
//print("cur:" + cur );
dir = cwd;

t = lengthOf(dir) == 0;
//print("t:" + t);

if (t)  {
   dir = cur;
}

//print("dir:" + dir );
args = "open=[" + dir + "] file=tif sort use";

//print( "\n" );
//print( "image format: " + format );
//print( "path: " + impath );
//print( "image: " + getDirectory("image") );
//print( "temp: " + getDirectory("temp") );
//print( "home: " + getDirectory("home") );
//print( "startup: " + getDirectory("startup") );
//print( "imagej: " + getDirectory("imagej") );
//print("current: " + getDirectory("current"));
//print("args: " + args)

// Auto contrast stack
run("Appearance...", "no auto ij menu=0 gui=1 16-bit=Automatic");

// Open image stack
run("Image Sequence...", args);

