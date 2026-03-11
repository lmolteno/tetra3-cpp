# Tetra3 C++
view the [tetra3 repo](https://github.com/esa/tetra3) for more info

## Intentions!
### Polar alignment
Main pain-point for setting up in the wops!
With plate-solving polar alignment should be a breeze.
This aspect will require a GPS to give accurate alt/az adjustments, but there may be another workaround. Could use the phone GNSS if we go for BLE integration, but I don't want any reliance on external hardware. Also - maybe plug in a ZWO to get super accurate solves through a scope? This adds a lot of complexity (i.e. how do we tune the exposure settings on the ZWO without visual feedback?), and I'm interested in the accuracy of this solver anyway

### ESP32 + super cheap camera star tracker
I did end up getting this to work - i could fit the southern sky onto the PSRAM. It could run solves fairly quickly and display the result on the OLED. The issue really came from the camera. The cheap OV cameras with parallel interfaces were very limited with how much control I had - and I didn't have great lenses (fairly soft).

### Pi Zero
Higher-cost but opens up to the MIPI cameras. 
