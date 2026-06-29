"""Shared colour utilities used by ete_profile.py."""


def color_gradient(start, end, mix=0.5):
    """Blend between start and end colors; return a CSS hex string.

    start / end can be:
      - a named string: 'white', 'black'
      - a hex string:   '#a1b2c3'
      - an RGB float tuple: (0.1, 0.5, 0.9)  — values in [0, 1]
    mix=0.0 → start, mix=1.0 → end.
    """
    def to_rgb(c):
        if isinstance(c, (list, tuple)):
            return tuple(float(x) for x in c[:3])
        named = {'white': (1., 1., 1.), 'black': (0., 0., 0.),
                 'red': (1., 0., 0.), 'blue': (0., 0., 1.)}
        if isinstance(c, str):
            if c.lower() in named:
                return named[c.lower()]
            h = c.lstrip('#')
            return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))
        return (1., 1., 1.)

    s = to_rgb(start)
    e = to_rgb(end)
    blended = tuple(s[i] + (e[i] - s[i]) * float(mix) for i in range(3))
    return '#{:02x}{:02x}{:02x}'.format(
        max(0, min(255, int(blended[0] * 255))),
        max(0, min(255, int(blended[1] * 255))),
        max(0, min(255, int(blended[2] * 255))),
    )
