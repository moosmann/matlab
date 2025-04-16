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

args = "open=[" + dir + "] file=tif sort";
//print("args: " + args)

// Auto contrast stack
run("Appearance...", "ij auto menu=0 gui=1 16-bit=Automatic");

// Open image stack
run("Image Sequence...", args);
