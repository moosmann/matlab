// Open image sequence in current directory from shell

dir = call("java.lang.System.getProperty", "user.dir");
cwd = getDirectory("current")
//print("cwd:" + cwd );
//print("dir:" + dir );

t = lengthOf(dir) == 0;
//print("t:" + t);

if (t)  {
   dir = cwd;
}
//print("dir:" + dir );

//args = "open=[" + dir + "] file=tif sort use";
list = getFileList(dir);
for (i=0; i<list.length; i++) {
    if (endsWith(list[i], "h5"))
      fn = list[i];
 }

//print( "\n" );
//print( "get file name: " + fn );
//print( "image format: " + format );
//print( "path: " + impath );
//print( "image: " + getDirectory("image") );
//print( "temp: " + getDirectory("temp") );
//print( "home: " + getDirectory("home") );
//print( "startup: " + getDirectory("startup") );
//print( "imagej: " + getDirectory("imagej") );
//print("current: " + getDirectory("current"));
args = "open=" + dir + "/" + fn + " use";
//print("args: " + args);

// Auto contrast stack
run("Appearance...", "no auto ij menu=0 gui=1 16-bit=Automatic");

// Open h5

run("HDF5...", args);
//run("Image Sequence...", args);



