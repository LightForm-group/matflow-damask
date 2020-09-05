from matflow.scripting import main_func
from damask_parse.utils import volume_element_from_2D_microstructure


@main_func
def volume_element_from_microstructure_image(microstructure_image, depth, image_axes):
    volume_element = volume_element_from_2D_microstructure(
        microstructure_image,
        depth,
        image_axes
    )
    return volume_element
