import os
import shutil
from nis import match

import pytest

from pasmopy import Text2Model
from pasmopy.construction.reaction_rules import DetectionError


def test_preprocessing():
    for i in ["1", "2"]:
        if os.path.isdir(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                f"typo_{i}",
            )
        ):
            shutil.rmtree(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    f"typo_{i}",
                )
            )


def test_typo_detection():
    for i in ["1", "2"]:
        typo = Text2Model(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                f"typo_{i}.txt",
            ),
        )
        with pytest.raises(DetectionError) as e:
            typo.convert()
        assert "Maybe: 'binds'." in str(e.value)


def test_cleanup():
    for i in ["1", "2"]:
        assert os.path.isdir(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                f"typo_{i}",
            )
        )
        shutil.rmtree(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                f"typo_{i}",
            )
        )
