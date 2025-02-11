
# Moon Phase Image Generator 🌕✨

This project is a Python-based Moon phase image generator that leverages NASA's SPICE toolkit (via SpiceyPy) to accurately compute the Moon’s ephemeris data, including its phase, illumination, and position relative to both the Earth and the Sun. The generated images include a white Moon with a realistic shadow overlay representing the terminator, all composed together to form the final Moon phase image.

---

## Features

- **Accurate Ephemeris Calculations**: Uses SPICE kernels to compute the Moon’s phase angle, illuminated fraction, distance, and more.
- **Realistic Illumination**: Computes a pixel-by-pixel shadow overlay based on the Moon’s local geometry.
- **Image Compositing**: Combines a clean white Moon image with a shadow layer for a visually appealing final result.
- **Enhanced Logging**: Clear, human-friendly logging messages with emojis to guide you through each step of the process.

---

## Prerequisites

Before running the script, make sure you have the following installed:

- **Python 3**  
- **Required Python Packages**:
  - `numpy`
  - `matplotlib`
  - `spiceypy`
  - `astropy`
  - `Pillow`

You can install the required packages using pip:

```bash
pip install numpy matplotlib spiceypy astropy pillow
```

---

## SPICE Kernels

This project requires a SPICE meta-kernel file named `moon.tm` in the same directory as the script. This meta-kernel should reference the necessary SPICE kernels (e.g., leapseconds, planetary ephemerides) to allow the computations.

---

## Usage

1. **Prepare the Environment**:
    - Ensure you have Python 3 installed.
    - Install the required Python packages.
    - Place the `moon.tm` file in the project directory.

2. **Customize (Optional)**:
    - **Observation Time**: Modify the `time_str` variable in the `main()` function to set your desired UTC observation time.
    - **Observer Location**: Change the observer location (e.g., `new_zealand` or `nl_location`) in the script if needed.
    - **Resolution**: Adjust the `resolution` variable to change the image output resolution.

3. **Run the Script**:

   Make the script executable and run it:

   ```bash
   chmod +x your_script.py
   ./your_script.py
   ```

   Or run it directly with Python:

    ```bash
    python your_script.py
    ```

4. **Output**:  
   The script will create an `output` directory containing:
    - `moon.png`: A white Moon image (a white circle with a black border) on a transparent background.
    - `shadow.png`: A shadow overlay image computed from the Moon’s illumination.
    - `final.png`: The final composite image with the Moon and its terminator (shadow overlay).

---

## Logging

The script utilizes Python's `logging` module with a custom, emoji-enhanced format to provide clear and fun output messages as it processes the data. This makes it easy to follow along with each step—from loading the SPICE kernels to generating and compositing the images.

---

## Customization

- **Observation Time**:  
  Change the `time_str` in the `main()` function.  
  Example:
  ```python
  time_str = '2026-05-31T18:30:00'
  ```

- **Observer Location**:  
  Modify the `location` variable in the `main()` function. The example uses New Zealand’s coordinates.  
  Example:
  ```python
  location = new_zealand
  ```

- **Image Resolution & Moon Radius**:  
  Adjust the `resolution` and `R` (Moon radius in arbitrary units) variables as desired for higher quality or differently scaled images.

---

## License

This project is open source. Feel free to modify, extend, and distribute it as needed.

---

## Acknowledgments

- **NASA SPICE Toolkit**: For providing precise planetary and spacecraft ephemeris data.
- **SpiceyPy**: The Python wrapper for the SPICE toolkit.
- **Astropy**: For time and coordinate transformations.
- **Matplotlib & Pillow**: For image creation and processing.

---

Enjoy generating stunning Moon phase images and may your projects always shine as bright as the Moon! 🌕🚀