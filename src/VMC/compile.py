"""Module containing the functions used to compile images into video format.

The main function of the module `compile_images()` is implemented within the `Protein.compile()` method.

"""

import os
import cv2 as cv
from tqdm import tqdm
from typing import List, Optional, Tuple, Union

def get_images(prefix: Optional[str],
               image_dir: str) -> Optional[List[str]]:

    """Returns the list of images to be compiled.

    Parameters
    ----------
    prefix : Optional[str]
        If none, all images in `image_dir` are returned.
        Otherwise, only images beginning with `prefix` are returned.

    image_dir : str
       The relative directory containing the images to be compiled.

    Returns
    -------
    Optional[List[str]]
        If no images are found in `image_dir` matching the (optional) `prefix`, None is returned.
        Otherwise, the list of image names to be compiled is returned.

    """
    # sort images in list by their scale number
    if prefix == None:
        images = sorted([image for image in os.listdir(image_dir)],
                        key=image_scale_number)

        if len(images) == 0:
            print(f"No images were found in {image_dir}, please try again")
            return None

    else:
        images = sorted([image for image in os.listdir(image_dir)
                         if image.startswith(prefix)],
                        key=image_scale_number)

        if len(images) == 0:
            print(f"No images were found in {image_dir} with 'prefix' {prefix}, please try again")
            return None

    return images

def image_scale_number(image_name: str) -> Union[int, str]:

    """Returns the scale number of the image from its name (if present).

    Used to sort images by their scale number when retrieving image names in `get_images()`.

    Parameters
    ----------
    image_name : str
        The name of the image.

    Returns
    -------
    Union[int, str]
        If the image contains a scale number, return the scale number as an integer.
        Otherwise, return the original image name as a string.

    """
    params = image_name.split("_")

    # tries to find scale number in image name
    try:
        scale_location = params.index("scale")
        if params[scale_location + 1].endswith(".png"):
            return int(params[scale_location + 1][:-4])
        else:
            return int(params[scale_location + 1])

    # if no scale number is found, returns image name
    except ValueError:
        return image_name

def get_dimensions(images: List[str], image_dir: str) -> Tuple[int, int]:

    """Returns the height and width of the image.

    Used to determine the dimensions of the final video.
    Only the first image in the list of images is passed to this function,
    as all images are assumed to have the same dimensions.

    Parameters
    ----------
    images : List[str]
        List of images to be compiled.

    image_dir : str
        The relative directory containing the images to be compiled.

    Returns
    -------
    Tuple[int, int]
        The height and width of the image.

    """
    frame = cv.imread(os.path.join(image_dir, images[0]))
    height, width = frame.shape[0], frame.shape[1]
    return height, width

def get_video_dir(video_name: Optional[str],
                  video_dir: str,
                  images: List[str]) -> str:

    """Returns the relative directory and name of the video to be saved.

    If `video_name` is specified, the specified name will be used.
    Otherwise, the video name is based off the name of the first image in the list of images.

    Parameters
    ----------
    video_name : Optional[str]
        Optional custom name for the video.

    video_dir : str
        Relative directory in which to save the video.

    images : List[str]
        List of images to be compiled.

    Returns
    -------
    str
        Relative directory and name of the video.

    """
    # if video name is specified, simply ensure that it ends in .mp4v
    if video_name != None:
        if not video_name.endswith(".mp4v"):
            video_name += ".mp4v"
        return f"{video_dir}/{video_name}"

    else:
        # if no video name is specified, base the video name off the first image name
        image_name = images[0]

        # tries to remove scale parameter from image name if present
        try:
            scale_location = image_name.index("scale")
            video_name = f"{image_name[:scale_location].rstrip('_')}.mp4v"
            return f"{video_dir}/{video_name}"

        # if no scale parameter is found, returns image name
        except ValueError:
            if image_name.endswith(".png"):
                video_name = f"{image_name[:-4]}.mp4v"
            else:
                video_name = f"{image_name}.mp4v"
            return f"{video_dir}/{video_name}"

def compile_images(prefix: Optional[str] = None,
                   video_name: Optional[str] = None,
                   fps: int = 2) -> None:

    """Compiles images into a .mp4v video file.

    The main function within this module, implemented within the `Protein.compile()` method.

    Parameters
    ----------
    prefix : Optional[str]
        If none, all images in `image_dir` are compiled.
        Otherwise, only images beginning with `prefix` are compiled.

    video_name : Optional[str]
        If specified, the video will be saved with this name.
        Otherwise, the video name is based on the name of the first image in the list of images.

    fps : int
        Frames per second of the video.
        A higher value corresponds to a faster playback speed of the video.

    """
    # specify relative image and video directories
    ##############################################
    image_dir = "./pymol_images"
    video_dir = "./pymol_videos"
    ##############################################

    # store all image files in a list
    images = get_images(prefix, image_dir)

    # checks that image list is not empty
    if images == None:
        return

    # get dimensions of photos
    height, width = get_dimensions(images, image_dir)

    # specify video name and location
    full_video_dir = get_video_dir(video_name, video_dir, images)

    # specify video options
    video = cv.VideoWriter(full_video_dir,
                           fourcc=0,
                           fps=fps,
                           frameSize=(width,height),
                           isColor=True)

    # write images to video
    for image_file in tqdm(images,
                           desc="Writing images to video",
                           unit="images"):

        image = cv.imread(os.path.join(image_dir, image_file))
        cv.putText(img=image,
                   text=image_file,
                   org=(100,2300),  # may have to be adjusted if cmd.png options are changed
                   fontFace=cv.FONT_HERSHEY_SIMPLEX,  # can be modified
                   fontScale=2,  # can be modified
                   color=(255,255,255),  # can be modified
                   thickness=3,  # can be modified
                   lineType=cv.LINE_AA,  # can be modified
                   bottomLeftOrigin=False)
        video.write(image)

    video.release()
    cv.destroyAllWindows()
    print(f"Animation '{full_video_dir}' has been generated successfully")

