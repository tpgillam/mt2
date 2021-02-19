from dataclasses import dataclass
from typing import Tuple, Union

import numpy


@dataclass(frozen=True)
class EllipseParams:
    """
    Represent an ellipse by the following form:

        c_xx x^2 + 2 c_xy x y + c_yy y^2 + 2 c_x x + 2 c_y y + c = 0

    Methods below follow the notation in:

        https://en.wikipedia.org/wiki/Matrix_representation_of_conic_sections
    """

    c_xx: float
    c_yy: float
    c_xy: float
    c_x: float
    c_y: float
    c: float

    @property
    def Aq(self) -> numpy.ndarray:
        """The matrix representation of the conic section."""
        return numpy.asarray(
            [
                [self.c_xx, self.c_xy, self.c_x],
                [self.c_xy, self.c_yy, self.c_y],
                [self.c_x, self.c_y, self.c],
            ]
        )

    @property
    def A33(self) -> numpy.ndarray:
        """The upper-left 2x2 section of Aq."""
        return self.Aq[:2, :2]

    @property
    def det_A33(self) -> float:
        """Determinant of A33."""
        return self.c_xx * self.c_yy - self.c_xy ** 2

    @property
    def det_Aq(self) -> float:
        """Determinant of Aq."""
        return numpy.linalg.det(self.Aq)

    @property
    def is_non_degenerate_real_ellipse(self) -> bool:
        """Returns true iff the equation represents a non-degenerate real ellipse."""
        if self.det_A33 < 0:
            return False
        if self.det_Aq == 0:
            return False
        factor = (self.c_xx + self.c_yy) * self.det_Aq
        return factor < 0

    @property
    def centre(self) -> Tuple[float, float]:
        """Get the co-ordinates of the centre of this conic section."""
        return (
            (self.c_xy * self.c_y - self.c_yy * self.c_x) / self.det_A33,
            (self.c_xy * self.c_x - self.c_xx * self.c_y) / self.det_A33,
        )

    def get_points(self, n_points: int) -> numpy.ndarray:
        """Return an array of shape (n_points, 2) of points on the ellipse."""
        centred_K = -self.det_Aq / self.det_A33
        transformed_A33 = self.A33 / centred_K
        L = numpy.linalg.cholesky(transformed_A33)

        theta = numpy.linspace(0, 2 * numpy.pi, n_points)
        # (n_points, 2)
        x_twiddle_prime = numpy.stack([numpy.cos(theta), numpy.sin(theta)], axis=1)
        x_twiddle = x_twiddle_prime @ numpy.linalg.inv(L.T).T
        return x_twiddle + numpy.asarray(self.centre)


def _make_ellipse_params(
    m_sq: float,  # Test parent squared mass
    mt_sq: float,  # Visible squared mass
    tx: float,  # Visible x-momentum
    ty: float,  # Visible y-momentum
    mq_sq: float,  # Invisible squared mass
    qx: float,  # Invisible x-momentum
    qy: float,  # Invisible y-momentum
) -> EllipseParams:
    """Construct an ellipse from physical quanties.

    This is just a translated version of `Lester::helper` from `lester_mt2_bisect_v7.h`
    """

    tx_sq = tx ** 2
    ty_sq = ty ** 2
    qx_sq = qx ** 2
    qy_sq = qy ** 2

    c_xx = +4.0 * mt_sq + 4.0 * ty_sq

    c_yy = +4.0 * mt_sq + 4.0 * tx_sq

    c_xy = -4.0 * tx * ty

    c_x = (
        -4.0 * mt_sq * qx
        - 2.0 * mq_sq * tx
        + 2.0 * m_sq * tx
        - 2.0 * mt_sq * tx
        + 4.0 * qy * tx * ty
        - 4.0 * qx * ty_sq
    )

    c_y = (
        -4.0 * mt_sq * qy
        - 4.0 * qy * tx_sq
        - 2.0 * mq_sq * ty
        + 2.0 * m_sq * ty
        - 2.0 * mt_sq * ty
        + 4.0 * qx * tx * ty
    )

    c = (
        -mq_sq * mq_sq
        + 2 * mq_sq * m_sq
        - m_sq * m_sq
        + 2 * mq_sq * mt_sq
        + 2 * m_sq * mt_sq
        - mt_sq * mt_sq
        + 4.0 * mt_sq * qx_sq
        + 4.0 * mt_sq * qy_sq
        + 4.0 * mq_sq * qx * tx
        - 4.0 * m_sq * qx * tx
        + 4.0 * mt_sq * qx * tx
        + 4.0 * mq_sq * tx_sq
        + 4.0 * qy_sq * tx_sq
        + 4.0 * mq_sq * qy * ty
        - 4.0 * m_sq * qy * ty
        + 4.0 * mt_sq * qy * ty
        - 8.0 * qx * qy * tx * ty
        + 4.0 * mq_sq * ty_sq
        + 4.0 * qx_sq * ty_sq
    )

    return EllipseParams(c_xx, c_yy, c_xy, c_x, c_y, c)


def make_ellipses(
    proposed_mt2: float,
    m_vis_1: Union[float, numpy.ndarray],
    px_vis_1: Union[float, numpy.ndarray],
    py_vis_1: Union[float, numpy.ndarray],
    m_vis_2: Union[float, numpy.ndarray],
    px_vis_2: Union[float, numpy.ndarray],
    py_vis_2: Union[float, numpy.ndarray],
    px_miss: Union[float, numpy.ndarray],
    py_miss: Union[float, numpy.ndarray],
    m_invis_1: Union[float, numpy.ndarray],
    m_invis_2: Union[float, numpy.ndarray],
) -> Tuple[EllipseParams, EllipseParams]:
    """
    Make a pair of ellipses in p1x, p1y; the intersection is the feasible region.
    """
    ellipse_1 = _make_ellipse_params(
        proposed_mt2 ** 2, m_vis_1 ** 2, -px_vis_1, -py_vis_1, m_invis_1 ** 2, 0, 0
    )
    ellipse_2 = _make_ellipse_params(
        proposed_mt2 ** 2,
        m_vis_2 ** 2,
        px_vis_2,
        py_vis_2,
        m_invis_2 ** 2,
        px_miss,
        py_miss,
    )
    return ellipse_1, ellipse_2
