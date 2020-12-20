/*
 * Generate bilinear registered images using TurboReg plugin to validate
 * the pystackreg registration/transformation implementation.
 */

var GOLDEN_RATIO = 0.5 * (sqrt(5.0) - 1.0);

var prefix = "pc12";
var path = File.directory();

var basename = prefix + "-unreg"

open(path + "/" + basename + ".tif");
n_images = nSlices();
run("Stack to Images");
selectWindow(basename + "-0001");
height = getHeight();
width = getWidth();

for (i=2; i<=n_images; i++) {
run("TurboReg",
	"-align "
	+ "-window " + basename + "-000" + i + " " // Source (window reference).
	+ "0 0 " + (width - 1) + " " + (height - 1) + " " // No cropping.
	+ "-window " + basename + "-0001" + " "// Target (window reference).
	+ "0 0 " + (width - 1) + " " + (height - 1) + " " // No cropping.
	+ "-bilinear "
	+ floor(width / 4 * GOLDEN_RATIO) + " " + floor(height / 4 * GOLDEN_RATIO) + " "
	+ floor(width / 4 * GOLDEN_RATIO) + " " + floor(height / 4 * GOLDEN_RATIO) + " "
	+ (floor(0.25 * GOLDEN_RATIO * width)) + " " + (height - (-floor(-0.25 * GOLDEN_RATIO * height))) + " "
	+ (floor(0.25 * GOLDEN_RATIO * width)) + " " + (height - (-floor(-0.25 * GOLDEN_RATIO * height))) + " "
	+ (width - (-floor(-0.25 * GOLDEN_RATIO * width))) + " " + (floor(0.25 * GOLDEN_RATIO * height)) + " "
	+ (width - (-floor(-0.25 * GOLDEN_RATIO * width))) + " " + (floor(0.25 * GOLDEN_RATIO * height)) + " "
	+ (width - (-floor(-0.25 * GOLDEN_RATIO * width))) + " " + (height - (-floor(-0.25 * GOLDEN_RATIO * height))) + " "
	+ (width - (-floor(-0.25 * GOLDEN_RATIO * width))) + " " + (height - (-floor(-0.25 * GOLDEN_RATIO * height))) + " "
	+ "-showOutput"); // In case -hideOutput is selected, the only way to

	selectWindow("Output");
	run("Stack to Images");
	selectWindow("Mask");
	rename("mask-reg-bilinear-0-" + (i-1));
	//close();

	selectWindow("Data");
	rename("reg-bilinear-0-" + (i-1));
	//saveAs("Tiff", path + prefix + "-reg-bilinear-0-" + (i-1) + ".tif");
	//close();
}

selectWindow(basename + "-0001");
rename("reg-bilinear-0-0");

for (i=2; i<=n_images; i++) {
	selectWindow(basename + "-000" + i);
	close();
}

run("Images to Stack", "name=Stack title=[mask-] use");
//saveAs("Tiff", path + prefix + "-reg-bilinear-mask.tif");

run("Images to Stack", "name=Stack title=[] use");
saveAs("Tiff", path + prefix + "-reg-bilinear.tif");
