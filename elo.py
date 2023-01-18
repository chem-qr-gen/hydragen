# implemented based on https://en.wikipedia.org/wiki/Elo_rating_system#Mathematical_details

def rate_single(rating_a: float, rating_b: float, score_a: float, score_b: float, max_adj: float = 32.0, ) -> (float, float):
    expected_score_a = 1/(1 + 10 ** ((rating_b - rating_a) / 400))
    expected_score_b = 1/(1 + 10 ** ((rating_a - rating_b) / 400))

    new_rating_a = rating_a + max_adj * (score_a - expected_score_a)
    new_rating_b = rating_b + max_adj * (score_b - expected_score_b)

    return (new_rating_a, new_rating_b)