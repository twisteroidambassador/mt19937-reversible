"""This is an educational-, research-, and CTF-oriented implementation of the
Mersenne Twister PRNG, providing features such as recovering internal state from
outputs, rewinding the internal state, or recovering the seed from a freshly
seeded internal state.

Copyright (C) 2021  twisteroid ambassador
Licensed under GPLv3
"""

import math
from typing import Optional, List


class MT19937_64:
    """Implements the 64-bit Mersenne Twister PRNG.
    Provides features such as recovering internal state from outputs, rewinding
    the internal state, or recovering the seed from a freshly seeded internal
    state.
    """

    """These are the algorithmic parameters for MT19937-64.
    Reference: https://en.wikipedia.org/wiki/Mersenne_Twister
    """
    # word size (bits)
    w = 64
    # degree of recurrence, also number of outputs required to clone
    n = 312
    m = 156
    r = 31
    a = 0xB5026F5AA96619E9
    u = 29
    d = 0x5555555555555555
    s = 17
    b = 0x71D67FFFEDA60000
    t = 37
    c = 0xFFF7EEE000000000
    l = 43
    f = 6364136223846793005

    def __init__(self, seed: Optional[int] = None):
        """Create an instance of the PRNG.
        If seed is not provided, the instance will have no internal state.
        Make sure to either seed, load state, or clone state from output before
        using.
        """

        """The internal state has length self.n. One random output is extracted
        from each state in order, and after extraction the state is changed by
        a twist operation. This is usually implemented in one of two ways: either
        do a mass twist on every state when all states have been used, or
        twist the just-used state after every extraction of random output.
        
        This implementation use the one-twist-per-output method, but with a
        twist (haha). A sane programmer would create a fixed states array,
        use a counter to point to the current state, twist the current state
        in place and increment the counter to move to the next state. Instead,
        we always keep the current state at the beginning of our states list,
        create the new twisted state at the end of the list and remove the
        current state from the beginning. This is obviously terrible for
        performance, but it also means we can use constant subscripts throughout
        the code, and makes it easier to reason about.
        """
        self._state: Optional[List[int]] = None
        if seed:
            self.seed(seed)

    def dump_state(self) -> List[int]:
        """Dump the internal state.
        The returned state can be used in load_state().
        """
        assert self._state is not None
        return self._state.copy()

    def load_state(self, state: List[int]) -> None:
        """Load the provided internal state into this instance."""
        assert len(state) == self.n
        self._state = state.copy()

    @staticmethod
    def _make_mask(bits: int) -> int:
        return (2 ** bits) - 1

    def _mask(self, num: int) -> int:
        """Keep the lower w bits."""
        return num & self._make_mask(self.w)

    def _upper_mask(self, num: int) -> int:
        """Keep the upper (w-r) bits of a w-bit number."""
        mask = self._make_mask(self.w - self.r) << self.r
        return num & mask

    def _lower_mask(self, num: int) -> int:
        """Keep the lower r bits."""
        return num & self._make_mask(self.r)

    def _times_A(self, x: int) -> int:
        """Simulate the multiply-by-A matrix operation."""
        xA = x >> 1
        if x & 1:
            xA = xA ^ self.a
        return xA

    def _invert_times_A(self, xA: int) -> int:
        """Undo the multiply-by-A matrix operation."""

        """Since x is shifted right one bit, its highest bit is always unset.
        Therefore by testing the highest bit we can know which branch in
        _times_A() has been taken, provided *a* has its highest bit set.
        Reference: https://www.ambionics.io/blog/php-mt-rand-prediction
        """
        assert self.a & (1 << (self.w - 1))
        if xA & (1 << (self.w - 1)):
            x = (xA ^ self.a) << 1 | 1
        else:
            x = xA << 1
        return x

    def _twist_one(self) -> None:
        """Twist the current state and move it to the end of the states list."""
        x = self._upper_mask(self._state[0]) | self._lower_mask(self._state[1])
        xA = self._times_A(x)
        new_state = self._state[self.m] ^ xA
        self._state.append(new_state)
        self._state.pop(0)

    def _untwist_one(self) -> None:
        """Move the last state to the beginning and untwist it."""
        xA = self._state[-1] ^ self._state[self.m - 1]
        x = self._invert_times_A(xA)
        last_state_high_bits = self._upper_mask(x)
        prev_xA = self._state[-2] ^ self._state[self.m - 2]
        prev_x = self._invert_times_A(prev_xA)
        last_state_low_bits = self._lower_mask(prev_x)
        last_state = last_state_low_bits | last_state_high_bits
        self._state.insert(0, last_state)
        self._state.pop(-1)

    def _temper(self, x: int) -> int:
        y1 = x ^ ((x >> self.u) & self.d)
        y2 = y1 ^ ((y1 << self.s) & self.b)
        y3 = y2 ^ ((y2 << self.t) & self.c)
        z = y3 ^ (y3 >> self.l)
        return z

    def _untemper(self, z: int) -> int:
        y3 = self._invert_xor_shift_and(z, self.l, self._make_mask(self.w))
        y2 = self._invert_xor_shift_and(y3, -self.t, self.c)
        y1 = self._invert_xor_shift_and(y2, -self.s, self.b)
        x = self._invert_xor_shift_and(y1, self.u, self.d)
        return x

    def _invert_xor_shift_and(self, num: int, shift_right: int, mask: int) -> int:
        """Undo the xor-with-shifted-and-masked-copy-of-self operation used in
        the tempering function.

        Provide a negative value for *shift_right* to shift left.

        Reference: https://occasionallycogent.com/inverting_the_mersenne_temper/index.html
        """
        rounds_required = math.ceil(self.w / abs(shift_right)) - 1
        intermediate = num
        for _ in range(rounds_required):
            if shift_right > 0:
                intermediate >>= shift_right
            else:
                intermediate <<= -shift_right
            intermediate = num ^ (intermediate & mask)
        return intermediate

    def seed(self, seed) -> None:
        """Seed the PRNG.
        *seed* must be a w-bit integer.
        """
        assert 0 <= seed < 2**self.w
        self._state = [seed]
        for i in range(1, self.n):
            self._state.append(self._mask(
                self.f * (self._state[-1] ^ (self._state[-1] >> (self.w - 2))) + i
            ))
        for _ in range(self.n):
            self._twist_one()

    def try_recover_seed(self) -> Optional[int]:
        """Recover the seed from a freshly-seeded state.

        If the PRNG has just been seeded, or its internal state is consistent
        with a just-seeded state, return the seed. Otherwise, return None.

        Reference: https://www.ambionics.io/blog/php-mt-rand-prediction
        """
        current_state = self.dump_state()
        try:
            self.rewind(self.n)
            """Even though all states are twisted once during seeding, simply
            untwisting all states will not put the seed back at _state[0].
            
            Conceptually, on an infinite state list:
            state[0] is the seed, state[1:n] are populated from the seed, and
            state[n:] are twisted from previous states:
            
            state[i+n] = twist(upper_bits(state[i]), lower_bits(state[i+1]), state[i+m])
            
            Therefore, state[0] only ever had its upper bits used in any
            twisting operation, and its lower bits cannot be recovered by
            un-twisting.
            
            On the other hands, state[1:n] can be recovered correctly, and it's
            easy to recover the seed from any one of these states.
            """
            for i in range(2, self.n):
                """Check whether the remaining states are consistent with being
                populated from the seed. If not, we're not in a freshly-seeded
                state."""
                if self._state[i] != self._mask(
                    self.f * (self._state[i-1] ^ (self._state[i-1] >> (self.w - 2))) + i
                ):
                    return None
            inv_f = pow(self.f, -1, 1 << self.w)
            u = self._mask(inv_f * (self._state[1] - 1))
            seed = u ^ (u >> (self.w - 2))
            assert self._upper_mask(seed) == self._upper_mask(self._state[0])
            return seed
        finally:
            self._state = current_state

    def get_next_random(self) -> int:
        """Generate one random output."""
        y = self._temper(self._state[0])
        self._twist_one()
        return y

    def rewind(self, rounds: int = 1) -> None:
        """Rewind the internal state for *rounds* iterations."""
        for _ in range(rounds):
            self._untwist_one()

    def clone_state_from_output_and_rewind(self, output: List[int]) -> None:
        """Clone internal states from a list of raw PRNG output.

        The list must have length of exactly self.n.

        After cloning, the internal state is at the beginning of the output list. In other words,
        self will begin to generate the output list, and the following statement will be True:
            output == [self.get_next_random() for _ in range(self.n)]
        """
        assert len(output) == self.n
        self._state = [self._untemper(z) for z in output]

    def clone_state_from_output(self, output: List[int]) -> None:
        """Clone internal states from a list of raw PRNG output.
        The list must have length of exactly self.n.
        After cloning, the internal state is at the end of the output list. In other words,
        if you take another Mersenne Twister PRNG, generate self.n outputs from it and
        feed into clone_states_from_output(), then the other PRNG and self will generate the
        same sequence of output from now on.
        """
        self.clone_state_from_output_and_rewind(output)
        for _ in range(self.n):
            self._twist_one()


class MT19937(MT19937_64):
    """Implements the 32-bit Mersenne Twister PRNG."""
    w = 32
    n = 624
    m = 397
    r = 31
    a = 0x9908B0DF
    u = 11
    d = 0xFFFFFFFF
    s = 7
    b = 0x9D2C5680
    t = 15
    c = 0xEFC60000
    l = 18
    f = 1812433253


def test_with_stdlib_random():
    import random
    stdlib_random = random.Random()
    mt = MT19937()
    output_for_cloning = [stdlib_random.getrandbits(mt.w) for _ in range(mt.n)]
    mt.clone_state_from_output(output_for_cloning)
    for _ in range(10000):
        assert stdlib_random.getrandbits(mt.w) == mt.get_next_random()
    print('Test successful')


if __name__ == '__main__':
    test_with_stdlib_random()
