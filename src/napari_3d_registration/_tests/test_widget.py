from napari_3d_registration import TimeRegistrationWidget, SpatialRegistrationWidget

# make_napari_viewer is a pytest fixture that returns a napari viewer object
# capsys is a pytest fixture that captures stdout and stderr output streams
def test_example_q_widget(make_napari_viewer, capsys):
    print("hello")