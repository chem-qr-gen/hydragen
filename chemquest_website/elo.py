# implemented based on https://en.wikipedia.org/wiki/Elo_rating_system#Mathematical_details

def rate_single(rating_a: float, rating_b: float, score_a: float, score_b: float, max_adj: float = 32.0) -> tuple[float, float]:
    expected_score_a = 1/(1 + 10 ** ((rating_b - rating_a) / 400))
    expected_score_b = 1/(1 + 10 ** ((rating_a - rating_b) / 400))

    new_rating_a = rating_a + max_adj * (score_a - expected_score_a)
    new_rating_b = rating_b + max_adj * (score_b - expected_score_b)

    return (new_rating_a, new_rating_b)

def calculate_new_elo(old_elo: float, question_elo: float, attempts: list[bool]) -> float:
    score = 0.0
    for attempt in attempts:
        if attempt == False:
            score -= 0.4
        else:
            score += 1.0
            break
    
    score = max(0.0, score)
    
    new_elo = rate_single(old_elo, question_elo, score, 1.0 - score)[0]
    return new_elo