"use strict";

var Filters = Filters || {};

////////////////////////////////////////////////////////////////////////////////
// General utility functions
////////////////////////////////////////////////////////////////////////////////

// Hardcoded Pi value
// const pi = 3.14159265359;
const pi = Math.PI;

// Constrain val to the range [min, max]
function clamp(val, min, max) {
  /* Shorthand for:
   * if (val < min) {
   *   return min;
   * } else if (val > max) {
   *   return max;
   * } else {
   *   return val;
   * }
   */
  return val < min ? min : val > max ? max : val;
}

// Extract vertex coordinates from a URL string
function stringToCoords(vertsString) {
  let centers = [];
  let coordStrings = vertsString.split("x");
  let coordsSoFar = 0;
  for (let i = 0; i < coordStrings.length; i++) {
    let coords = coordStrings[i].split("y");
    let x = parseInt(coords[0]);
    let y = parseInt(coords[1]);
    if (!isNaN(x) && !isNaN(y)) {
      centers.push({ x: x, y: y });
    }
  }

  return centers;
}

// Blend scalar start with scalar end. Note that for image blending,
// end would be the upper layer, and start would be the background
function blend(start, end, alpha) {
  return start * (1 - alpha) + end * alpha;
}

// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 73 lines of code.
// ----------- STUDENT CODE END ------------

////////////////////////////////////////////////////////////////////////////////
// Filters
////////////////////////////////////////////////////////////////////////////////

// You've already implemented this in A0! Feel free to copy your code into here
Filters.fillFilter = function (image, color) {
  for (let x = 0; x < image.width; ++x) {
    for (let y = 0; y < image.height; ++y) {
      image.setPixel(x, y, color);
    }
  }
  return image;
};

// You've already implemented this in A0! Feel free to copy your code into here
Filters.brushFilter = function (image, radius, color, vertsString) {
  // centers is an array of (x, y) coordinates that each defines a circle center

  const centers = stringToCoords(vertsString);

  // draw a filled circle centered at every location in centers[].
  // radius and color are specified in function arguments.
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 10 lines of code.

  for (let cord in centers) {
    for (let x = 0; x < image.width; ++x) {
      for (let y = 0; y < image.height; ++y) {
        if (
          Math.sqrt(
            Math.abs(centers[cord].x - x) ** 2 +
              Math.abs(centers[cord].y - y) ** 2
          ) < radius
        ) {
          /*   if (
          Math.abs(centers[cord].x - x) + Math.abs(centers[cord].y - y) <
          radius
        ) {*/
          image.setPixel(x, y, color);
        }
      }
    }
  }
  return image;
};

// You've already implemented this in A0! Feel free to copy your code into here
/* Filters.softBrushFilter = function (
  image,
  radius,
  color,
  alpha_at_center,
  vertsString
) {
  // centers is an array of (x, y) coordinates that each defines a circle center
  const centers = stringToCoords(vertsString);
  // console.log(centers);

  // draw a filled circle with opacity equals to alpha_at_center at the center of each circle
  // the opacity decreases linearly along the radius and becomes zero at the edge of the circle
  // radius and color are specified in function arguments.
  // ----------- STUDENT CODE BEGIN ------------
  for (let i in centers) {
    let el = centers[i];
    for (let x = 0; x < image.width; ++x) {
      for (let y = 0; y < image.height; ++y) {
        // let beta = alpha_at_center;
        // if (beta < 0) beta = 0;
        let orig = image.getPixel(x, y);
        //(x,y)ի հեռավորությունը (el.x, el.y)

        let heravorutyun = Math.sqrt((x - el.x) ** 2 + (y - el.y) ** 2);
        let haraberakan = heravorutyun / radius;

        if (haraberakan < 1) {
          let alpha = (1 - haraberakan) * alpha_at_center;
          image.setPixel(
            x,
            y,
            orig.multipliedBy(1 - alpha).plus(color.multipliedBy(alpha))
          );
        }
        //orig: նկարի միջի նախկին գույն
        //color: փափուկ վրձնի գույնր
        //orig*(1-beta)+color*beta;

        //orig: [1,0.9,1]
        //brush: [0,0,0]
        //միջին։ [(1+0)/2, (0.9+0)/2, (1+0)/2]
      }
    }

    // ----------- Our reference solution uses 20 lines of code.
    // ----------- STUDENT CODE END ------------

    return image;
  }
};  */

Filters.softBrushFilter = function( image, radius, color, alpha_at_center, vertsString ) {
    // centers is an array of (x, y) coordinates that each defines a circle center
    const centers = stringToCoords(vertsString);
    for (let i in centers) {
    let el = centers[i];
    for (let x = 0; x < image.width; x++) {
      for (let y = 0; y < image.height; y++) {
        // let beta = alpha_at_center;
        // if (beta < 0) beta = 0;
        let orig = image.getPixel(x, y);
        // (x,y)ի հեռավորությունը (el.x, el.y)
        let heravorutyun = Math.sqrt((x-el.x) ** 2 + (y-el.y) ** 2);
        let haraberakan = heravorutyun/radius;

        if (haraberakan < 1) {
            let alpha = (1 - haraberakan) * alpha_at_center;
            image.setPixel(x, y, orig.multipliedBy(1-alpha).plus(color.multipliedBy(alpha)));
        }

        // orig: [1,0.9,1]
        // brush: [0,0, 0]
        // միջին։  [(1+0)/2, (0.9+0)/2, (1+0)/2]
        //        [1*0.9 + 0*0.1, 0.9*0.9+0*0.1, 1*0.9+0*0.1]
        //        [1*beta + 0*(1-beta), 0.9*beta+0*(1-beta), 1*beta+0*(1-beta)]

        // orig: նկարի միջի նախկին գույնը
        // color: փափուկ վրձինի գույնը
        // orig * (1-beta)+color*beta
      }
    }
    }
    // ----------- STUDENT CODE END ------------

    return image;
};

// Ratio is a value in the domain [-1, 1]. When ratio is < 0, linearly blend the image
// with black. When ratio is > 0, linearly blend the image with white. At the extremes
// of -1 and 1, the image should be completely black and completely white, respectively.
Filters.brightnessFilter = function (image, ratio) {
  let alpha, dirLuminance;
  if (ratio < 0.0) {
    alpha = Math.abs(ratio);
    dirLuminance = 0; // blend with black
  } else {
    alpha = 1 - Math.abs(ratio);
    dirLuminance = 1; // blend with white
  }

  for (var x = 0; x < image.width; x++) {
    for (var y = 0; y < image.height; y++) {
      var pixel = image.getPixel(x, y);

      pixel.data[0] = alpha * pixel.data[0] + (1 - alpha) * dirLuminance;
      pixel.data[1] = alpha * pixel.data[1] + (1 - alpha) * dirLuminance;
      pixel.data[2] = alpha * pixel.data[2] + (1 - alpha) * dirLuminance;

      image.setPixel(x, y, pixel);
    }
  }

  return image;
};

// Reference at this:
//      https://en.wikipedia.org/wiki/Image_editing#Contrast_change_and_brightening
// value = (value - 0.5) * (tan ((contrast + 1) * PI/4) ) + 0.5;
// Note that ratio is in the domain [-1, 1]



/* Filters.contrastFilter = function (image, ratio) {
  if (ratio < 0.0) ratio = ratio * (1.0 + ratio);
  else ratio = ratio + (1 - ratio) * ratio;

  for (var x = 0; x < image.width; x++) {
    for (var y = 0; y < image.height; y++) {
      let pixel = image.getPixel(x, y);
      pixel.data[0] =
        (pixel.data[0] - 0.5) * Math.tan(((ratio + 1) * Math.PI) / 4) + 0.5;
      pixel.data[1] =
        (pixel.data[1] - 0.5) * Math.tan(((ratio + 1) * Math.PI) / 4) + 0.5;
      pixel.data[2] =
        (pixel.data[2] - 0.5) * Math.tan(((ratio + 1) * Math.PI) / 4) + 0.5;

      image.setPixel(x, y, pixel);
    }
  }
  //Gui.alertOnce("contrastFilter is not implemented yet");
  return image;
};

// Note that the argument here is log(gamma)
Filters.gammaFilter = function (image, logOfGamma) {
  const gamma = Math.exp(logOfGamma);

  for (let x = 0; x < image.width; x++) {
    for (let y = 0; y < image.height; y++) {
      let pixel = image.getPixel(x, y);
      pixel.data[0] = pixel.data[0] * gamma;
      pixel.data[1] = pixel.data[1] * gamma;
      pixel.data[2] = pixel.data[2] * gamma;
      image.setPixel(x, y, pixel);
    }
  }
  //Gui.alertOnce("gammaFilter is not implemented yet");
  return image;
};  */



Filters.contrastFilter = function( image, ratio ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 14 lines of code.
    let radius = (ratio+1) * 15;
    // radius is in range  [0, 30]
    let copy = image.copyImg();
    for (var x = 0; x < image.width; x++) {
        for (var y = 0; y < image.height; y++) {
            var pixel = image.getPixel(x, y);
            var sum = new Pixel(0,0,0);

            // sum.plus(image.getPixel(x-1, y-1));
            // sum.plus(image.getPixel(x-1, y));
            // sum.plus(image.getPixel(x-1, y+1));
            // sum.plus(image.getPixel(x, y-1));
            // sum.plus(image.getPixel(x, y));
            // sum.plus(image.getPixel(x, y+1));
            // sum.plus(image.getPixel(x+1, y-1));
            // sum.plus(image.getPixel(x+1, y));
            // sum.plus(image.getPixel(x+1, y+1));

            // [0,255,0,0,255]
            // [255/3, 255/3, 255/3, 255/3]

            // [255/3,255,0,0,255]
            //radius
           
           
             /*let filter = [[1,1,1], 
                          [1,1,1],
                          [1,1,1]]; */

          /*  let filter = [[-1,0,-1], 
                          [-1,6,-1],
                          [-1,0,-1]];
*/
                          
            /*let filter =     [[-1,-1,-1], 
                          [-1,8,-1],
                          [-1,-1,-1]];*/
                          
         /*  let filter = [[-1.5,1.5,0], 
                           [-1.5,1.5,0],
                           [-1.5,1.5,0]];            */
             let filter = [[0,1,0], 
                           [0,0,0], 
                           [0,0,0]]; 
            radius = 1;
            for (let i = x-radius; i <= x+radius; i++) {
                for (let j = y-radius; j<= y+radius; j++) {
                    // console.log(filter, filter[i+radius-x], filter[i+radius-x][j+radius-y])
                    sum = sum.plus(copy.getPixel(i,j).multipliedBy(filter[j+radius-y][i+radius-x] ));
                }
            }
            // 1 , (1 + 2*1) ** 2 = 3 ** 2 = 9
           // sum = sum.dividedBy((1 + 2*radius) ** 2);

            image.setPixel(x, y, sum);
        }
    }
    // ----------- STUDENT CODE END ------------
    return image;

};




/*
 * The image should be perfectly clear up to innerRadius, perfectly dark
 * (black) at outerRadius and beyond, and smoothly increase darkness in the
 * circular ring in between. Both are specified as multiples of half the length
 * of the image diagonal (so 1.0 is the distance from the image center to the
 * corner).
 *
 * Note that the vignette should still form a perfect circle!
 */
Filters.vignetteFilter = function (image, innerR, outerR) {
  // Let's ensure that innerR is at least 0.1 smaller than outerR
  innerR = clamp(innerR, 0, outerR - 0.1);
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 17 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("vignetteFilter is not implemented yet");
  return image;
};

/*
 * You will want to build a normalized CDF of the L channel in the image.
 */
Filters.histogramEqualizationFilter = function (image) {
  
  //Gui.alertOnce("histogramEqualizationFilter is not implemented yet");
  return image;
};



/* false
*
*  Filters.histogramEqualizationFilter = function(img, x1, y1, x2, y2, num_bins) {
*    if( num_bins == undefined )
*        num_bins = 256;
*    var h = img.h, w = img.w;
*    var hist = [];
*    var i, x, y, idx, val;
*    // initialize the histogram
*    for(i=0;i<num_bins;++i)
*        hist[i] = 0;
*    // loop over every single pixel
*    for(y=y1,idx=0;y<y2;++y) {
*        for(x=x1;x<x2;++x,idx+=4) {
*            // figure out which bin it is in
*            val = Math.floor((img.data[idx] / 255.0) * (num_bins-1));
*            ++hist[val];
*        }
*   }
    return hist;
}
*/




// Set each pixel in the image to its luminance
Filters.grayscaleFilter = function (image) {
  for (let x = 0; x < image.width; x++) {
    for (let y = 0; y < image.height; y++) {
      const pixel = image.getPixel(x, y);
      const luminance =
        0.2126 * pixel.data[0] +
        0.7152 * pixel.data[1] +
        0.0722 * pixel.data[2];
      pixel.data[0] = luminance;
      pixel.data[1] = luminance;
      pixel.data[2] = luminance;

      image.setPixel(x, y, pixel);
    }
  }

  return image;
};

// Adjust each channel in each pixel by a fraction of its distance from the average
// value of the pixel (luminance).
// See: http://www.graficaobscura.com/interp/index.html
Filters.saturationFilter = function (image, ratio) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 13 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("saturationFilter is not implemented yet");
  return image;
};

// Apply the Von Kries method: convert the image from RGB to LMS, divide by
// the LMS coordinates of the white point color, and convert back to RGB.
Filters.whiteBalanceFilter = function (image, white) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 23 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("whiteBalanceFilter is not implemented yet");
  return image;
};

// This is similar to the histogram filter, except here you should take the
// the CDF of the L channel in one image and
// map it to another
//
Filters.histogramMatchFilter = function (image, refImg) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 58 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("histogramMatchFilter is not implemented yet");
  return image;
};

// Convolve the image with a gaussian filter
Filters.gaussianFilter = function (image, sigma) {
  // note: this function needs to work in a new copy of the image
  //       to avoid overwriting original pixels values needed later
  // create a new image with the same size as the input image
  var newImg = image.createImg(image.width, image.height);
  // the filter window will be [-winR, winR] for a total diameter of roughly Math.round(3*sigma)*2+1;
  const winR = Math.round(sigma * 3);
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 54 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("gaussianFilter is not implemented yet");
  return newImg;
};

/*
 * First the image with the edge kernel and then add the result back onto the
 * original image.
 */
Filters.sharpenFilter = function (image) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 29 lines of code.
  // ----------- STUDENT CODE END ----------
  let radius = 1;
    // radius is in range  [0, 30]
    let copy = image.copyImg();
    for (var x = 0; x < image.width; x++) {
        for (var y = 0; y < image.height; y++) {
            var pixel = image.getPixel(x, y);
            var sum = new Pixel(0,0,0);
                 
             let filter = [[-1,-1,-1], 
                           [-1,9,-1], 
                           [-1,-1,-1]]; 
            radius = 1;
            for (let i = x-radius; i <= x+radius; i++) {
                for (let j = y-radius; j<= y+radius; j++) {
                    // console.log(filter, filter[i+radius-x], filter[i+radius-x][j+radius-y])
                    sum = sum.plus(copy.getPixel(i,j).multipliedBy(filter[j+radius-y][i+radius-x] ));
                }
            }
            // 1 , (1 + 2*1) ** 2 = 3 ** 2 = 9
           // sum = sum.dividedBy((1 + 2*radius) ** 2);

            image.setPixel(x, y, new Pixel(0,0,0).minus(sum));
        }
    }


  //Gui.alertOnce("sharpenFilter is not implemented yet");
  return image;
};

/*
 * Convolve the image with the edge kernel from class. You might want to define
 * a convolution utility that convolves an image with some arbitrary input kernel
 *
 * For this filter, we recommend inverting pixel values to enhance edge visualization
 */
Filters.edgeFilter = function (image) {
  image = Filters.grayscaleFilter(image);
  let radius = 1;
    // radius is in range  [0, 30]
    let copy = image.copyImg();
    for (var x = 0; x < image.width; x++) {
        for (var y = 0; y < image.height; y++) {
            var pixel = image.getPixel(x, y);
            var sum = new Pixel(0,0,0);

             let filter = [[-1,-1,-1], 
                           [-1,8,-1], 
                           [-1,-1,-1]]; 
            radius = 1;
            for (let i = x-radius; i <= x+radius; i++) {
                for (let j = y-radius; j<= y+radius; j++) {
                    // console.log(filter, filter[i+radius-x], filter[i+radius-x][j+radius-y])
                    sum = sum.plus(copy.getPixel(i,j).multipliedBy(filter[j+radius-y][i+radius-x] ));
                }
            }
            // 1 , (1 + 2*1) ** 2 = 3 ** 2 = 9
           // sum = sum.dividedBy((1 + 2*radius) ** 2);

            image.setPixel(x, y, new Pixel(1,1,1).minus(sum));
        }
    }
    // ----------- STUDENT CODE END ------------
    return image;

};

 
// Set a pixel to the median value in its local neighbor hood. You might want to
// apply this seperately to each channel.
Filters.medianFilter = function (image, winR) {
  // winR: the window will be  [-winR, winR];
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 34 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("medianFilter is not implemented yet");
  return image;
};

// Apply a bilateral filter to the image. You will likely want to reference
// precept slides, lecture slides, and the assignments/examples page for help.
Filters.bilateralFilter = function (image, sigmaR, sigmaS) {
  // reference: https://en.wikipedia.org/wiki/Bilateral_filter
  // we first compute window size and preprocess sigmaR
  const winR = Math.round((sigmaR + sigmaS) * 1.5);
  sigmaR = sigmaR * (Math.sqrt(2) * winR);

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 49 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("bilateralFilter is not implemented yet");
  return image;
};

// Conver the image to binary
Filters.quantizeFilter = function (image) {
  // convert to grayscale
  image = Filters.grayscaleFilter(image);

  // use center color
  for (let i = 0; i < image.height; i++) {
    for (let j = 0; j < image.width; j++) {
      const pixel = image.getPixel(j, i);
      let rand = Math.random();
      for (let c = 0; c < 3; c++) {
        pixel.data[c] = Math.round(pixel.data[c] + rand);
      }
      pixel.clamp();
      image.setPixel(j, i, pixel);
    }
  }
  return image;
};

// To apply random dithering, first convert the image to grayscale, then apply
// random noise, and finally quantize
Filters.randomFilter = function (image) {
  // convert to grayscale
  image = Filters.grayscaleFilter(image);


  for (let i = 0; i < image.height; i++) {
    for (let j = 0; j < image.width; j++) {
let rand = Math.random() - 0.5;
const pixel = image.getPixel(j,i);
for(let c = 0; c < 3; c++){
  pixel.data[c] = Math.random(pixel.data[c] + rand);
}
pixel.clamp();
image.setPixel(j,i,pixel);
    }}

  return image;
};

// Apply the Floyd-Steinberg dither with error diffusion
Filters.floydFilter = function (image) {
  // convert to grayscale
 image = Filters.grayscaleFilter(image);

  for(let i = 0; i < image.height; i++){
    for(let j = 0; j < image.width; j++){
      const pixel = image.getPixel(j,i);
      let new_val = Math.round(pixel.data[0]);

      let error = new Pixel(0,0,0);
      for(let k = 0; k < 3; k++){
        error.data[k] = pixel.data[0]-new_val;
      }

      for(let c = 0; c< 3; c++){
        pixel.data[c] = new_val;
      }

      pixel.clamp();
      image.setPixel(j,i,pixel);

      
      if(i == 0 && j ==0 ) console.log(error);
      image.setPixel(j+1,i,image.getPixel(j+1,i).plus(error.multipliedBy(7/16)));
      image.setPixel(j-1,i+1,image.getPixel(j-1,i+1).plus(error.multipliedBy(5/16)));
      image.setPixel(j,i+1,image.getPixel(j,i+1).plus(error.multipliedBy(3/16)));
      image.setPixel(j+1,i+1,image.getPixel(j+1,i+1).plus(error.multipliedBy(1/16)));

    }
  }
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 24 lines of code.
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("floydFilter is not implemented yet");
  return image;
};

// Apply ordered dithering to the image. We recommend using the pattern from the
// examples page and precept slides.
Filters.orderedFilter = function (image) {
  // convert to gray scale
  image = Filters.grayscaleFilter(image);

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 30 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("orderedFilter is not implemented yet");
  return image;
};

// Implement bilinear and Gaussian sampling (in addition to the basic point sampling).
// This operation doesn't appear on GUI and should be used as a utility function.
// Call this function from filters that require sampling (e.g. scale, rotate)
Filters.samplePixel = function (image, x, y, mode) {
  if (mode == "bilinear") {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 19 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("bilinear sampling is not implemented yet");
  } else if (mode == "gaussian") {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 38 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("gaussian sampling is not implemented yet");
  } else {
    // point sampling
    y = Math.max(0, Math.min(Math.round(y), image.height - 1));
    x = Math.max(0, Math.min(Math.round(x), image.width - 1));
    return image.getPixel(x, y);
  }
};

// Translate the image by some x, y and using a requested method of sampling/resampling
Filters.translateFilter = function (image, x, y, sampleMode) {
  // Note: set pixels outside the image to RGBA(0,0,0,0)
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 21 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("translateFilter is not implemented yet");
  return image;
};

// Scale the image by some ratio and using a requested method of sampling/resampling
Filters.scaleFilter = function (image, ratio, sampleMode) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 19 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("scaleFilter is not implemented yet");
  return image;
};

// Rotate the image by some angle and using a requested method of sampling/resampling
Filters.rotateFilter = function (image, radians, sampleMode) {
  // Note: set pixels outside the image to RGBA(0,0,0,0)
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 30 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("rotateFilter is not implemented yet");
  return image;
};

// Swirl the filter about its center. The rotation of the swirl should be in linear increase
// along the radial axis up to radians
Filters.swirlFilter = function (image, radians, sampleMode) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 27 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("swirlFilter is not implemented yet");
  return image;
};

// Set alpha from luminance
Filters.getAlphaFilter = function (backgroundImg, foregroundImg) {
  for (let i = 0; i < backgroundImg.height; i++) {
    for (let j = 0; j < backgroundImg.width; j++) {
      const pixelBg = backgroundImg.getPixel(j, i);
      const pixelFg = foregroundImg.getPixel(j, i);
      const luminance =
        0.2126 * pixelFg.data[0] +
        0.7152 * pixelFg.data[1] +
        0.0722 * pixelFg.data[2];
      pixelBg.a = luminance;
      backgroundImg.setPixel(j, i, pixelBg);
    }
  }

  return backgroundImg;
};

// Composites the foreground image over the background image, using the alpha
// channel of the foreground image to blend two images.
Filters.compositeFilter = function (backgroundImg, foregroundImg) {
  // Assume the input images are of the same sizes.
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 14 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("compositeFilter is not implemented yet");
  return backgroundImg;
};

// Morph two images according to a set of correspondance lines
Filters.morphFilter = function( initialImg, finalImg, alpha, sampleMode, linesFile ) {
  const lines = Parser.parseJson( "images/" + linesFile );

  // The provided linesFile represents lines in a flipped x, y coordinate system
  //  (i.e. x for vertical direction, y for horizontal direction).
  // Therefore we first fix the flipped x, y coordinates here.
  for (let i = 0; i < lines.initial.length; i++) {
      [lines.initial[i].x0, lines.initial[i].y0] = [lines.initial[i].y0, lines.initial[i].x0];
      [lines.initial[i].x1, lines.initial[i].y1] = [lines.initial[i].y1, lines.initial[i].x1];
      [lines.final[i].x0, lines.final[i].y0] = [lines.final[i].y0, lines.final[i].x0];
      [lines.final[i].x1, lines.final[i].y1] = [lines.final[i].y1, lines.final[i].x1];
  }

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 86 lines of code.
  let p = 0.5;
  let a = 0.01;
  let b = 2;
  let alphaImg = initialImg.createImg(initialImg.width, initialImg.height);
  let alphaLines = [];
  for (let i = 0; i < lines.initial.length; i++) {
      alphaLines[i] = {};
      alphaLines[i].x0 = (1-alpha) * lines.initial[i].x0 + alpha * lines.final[i].x0;
      alphaLines[i].y0 = (1-alpha) * lines.initial[i].y0 + alpha * lines.final[i].y0;
      alphaLines[i].y1 = (1-alpha) * lines.initial[i].y1 + alpha * lines.final[i].y1;
      alphaLines[i].x1 = (1-alpha) * lines.initial[i].x1 + alpha * lines.final[i].x1;
  }

  //                          dest_lines, src_lines
  function actuallyMorph(img, lines0, lines1){
      /*
      For each pixel X in the destination
      DSUM = (0,0)
      weightsum = 0
      For each line Pi Qi
          calculate u,v based on Pi Qi
          calculate X'i based on u,v       and Pi'Qi'
          calculate displacement Di = Xi' - Xi for this line
          dist = shortest distance from X to Pi Qi
          weight = (lengthp / (a + dist))b
          DSUM += Di *    weight
          weightsum += weight
      X' = X + DSUM / weightsum
      destinationImage(X) = sourceImage(X')
      */
      let morphImg = img.createImg(img.width, img.height);
      for (let y = 0; y < img.height; y++) {
          for (let x = 0 ;x < img.width; x++) {
              // DSUM=(0,0)
              let sumX = 0;
              let sumY = 0;
              let weightsum = 0;
              for (let l = 0; l < lines0.length; l++) {
                  // (X-P) * (Q-P) = (x-Px) * (Qx-Px) + (y-Py)*(Qy-Py) = (x-x0) * (x1-x0) + (y-y0)*(y1-y0) // all points in l0
                  let V = lines0[l];// v as in vector as in point of 2 coords
                  let Vp = lines1[l];
                  // u = (X-V)
                  let u = ( (x-V.x0) * (V.x1-V.x0) + (y-V.y0)*(V.y1-V.y0) )/V.lenSq;
                  // perp(x,y) = (-y,x)
                  let v = ( (x-V.x0) * -1 * (V.y1-V.y0) + (y-V.y0)*(V.x1-V.x0) )/Math.sqrt(V.lenSq);

                  let xp = Vp.x0 + u * (Vp.x1-Vp.x0) + v *  -1 * (Vp.y1-Vp.y0) / Math.sqrt(Vp.lenSq);
                  let yp = Vp.y0 + u * (Vp.y1-Vp.y0) + v * (Vp.x1-Vp.x0) / Math.sqrt(Vp.lenSq);


                  // v = shortest distance from X to PQ
                  let dist = Math.abs(v)
                  if(u < 0)
                      dist = Math.sqrt((V.x0 - x) ** 2 + (V.y0 - y) ** 2);
                  else if (u > 1)
                      dist = Math.sqrt((V.x1 - x) ** 2 + (V.y1 - y) ** 2);

                  let w = (V.lenSq ** (p/2))/(a+dist);
                  w = w ** b;

                  // Di = Xi' - Xi
                  let dx = xp - x;
                  let dy = yp - y;
                  sumX += dx * w;
                  sumY += dy * w;
                  weightsum += w;

              }
          let xp = x + sumX/weightsum;
          let yp = y + sumY/weightsum;

          
          let pixel = Filters.samplePixel(img, xp, yp, sampleMode);
          morphImg.setPixel(x,y, pixel);
          }
      }

      return morphImg;
  };
  lines.initial.map(o => o.lenSq = (o.x0-o.x1) ** 2 + (o.y0-o.y1) ** 2);
  lines.final.map(o => o.lenSq = (o.x0-o.x1) ** 2 + (o.y0-o.y1) ** 2);
  alphaLines.map(o => o.lenSq = (o.x0-o.x1) ** 2 + (o.y0-o.y1) ** 2);

  let initMorph = actuallyMorph(initialImg, alphaLines, lines.initial);
  let finiMorph = actuallyMorph(finalImg, alphaLines, lines.final);


  for (let i = 0; i < initialImg.height; i++) {
      for (let j = 0; j < initialImg.width; j++) {
          const pixelInit = initMorph.getPixel(j, i);
          const pixelFini = finiMorph.getPixel(j, i);
          const pixel = pixelInit.multipliedBy(1-alpha).plus(pixelFini.multipliedBy(alpha))
          alphaImg.setPixel(j, i, pixel);
      }
  }
  // ----------- STUDENT CODE END ------------
  return alphaImg;
};

// Use k-means to extract a pallete from an image
Filters.paletteFilter = function (image, colorNum) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 82 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("paletteFilter is not implemented yet");
  return image;
};

// Read the following paper and implement your own "painter":
//      http://mrl.nyu.edu/publications/painterly98/hertzmann-siggraph98.pdf
Filters.paintFilter = function (image, value) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 54 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("paintFilter is not implemented yet");
  return image;
};

/*
 * Read this paper for background on eXtended Difference-of-Gaussians:
 *      http://www.cs.princeton.edu/courses/archive/spring19/cos426/papers/Winnemoeller12.pdf
 * Read this paper for an approach that develops a flow field based on a bilateral filter
 *      http://www.cs.princeton.edu/courses/archive/spring19/cos426/papers/Kang09.pdf
 */
Filters.xDoGFilter = function (image, value) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 61 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("xDoGFilter is not implemented yet");
  return image;
};

// You can use this filter to do whatever you want, for example
// trying out some new idea or implementing something for the
// art contest.
// Currently the 'value' argument will be 1 or whatever else you set
// it to in the URL. You could use this value to switch between
// a bunch of different versions of your code if you want to
// code up a bunch of different things for the art contest.
Filters.customFilter = function (image, value) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 0 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("customFilter is not implemented yet");
  return image;
};
