// Open image sequence in current directory from shell

cwd = call("java.lang.System.getProperty", "user.dir");
cur = getDirectory("current")
//print("cwd:" + cwd );
//print("cur:" + cur );
dir = cur;

t = lengthOf(dir) == 0;
//print("t:" + t);

if (t)  {
   dir = cwd;
}

//print("dir:" + dir );
args = "open=[" + dir + "] file=tif sort";

// Auto contrast stack
run("Appearance...", "auto ij menu=0 gui=1 16-bit=Automatic");

// Open image stack
run("Image Sequence...", args);

