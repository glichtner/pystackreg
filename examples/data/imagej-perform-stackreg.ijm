/*
 * Generate registered images using StackReg plugin to validate
 * the pystackreg registration/transformation implementation.
 */

var plugin_name = "StackReg";
var prefix = "pc12";
var path = File.directory();


open(path + "/" + prefix + "-unreg.tif");
run(plugin_name, "transformation=Translation");
saveAs("Tiff", path + "/" + prefix + "-reg-translation.tif");
close();

open(path + "/" + prefix + "-unreg.tif");
run(plugin_name, "transformation=[Rigid Body]");
saveAs("Tiff", path + "/" + prefix + "-reg-rigid-body.tif");
close();

open(path + "/" + prefix + "-unreg.tif");
run(plugin_name, "transformation=[Scaled Rotation]");
saveAs("Tiff", path + "/" + prefix + "-reg-scaled-rotation.tif");
close();

open(path + "/" + prefix + "-unreg.tif");
run(plugin_name, "transformation=[Affine]");
saveAs("Tiff", path + "/" + prefix + "-reg-affine.tif");
close();
